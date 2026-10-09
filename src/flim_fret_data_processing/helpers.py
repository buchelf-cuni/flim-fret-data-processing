"""Shared functions for the FLIM phasor pipeline (notebooks 01-05)."""
import contextlib
import os
import re
import struct
import zlib
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from scipy import ndimage
from scipy.ndimage import gaussian_filter1d
from scipy.optimize import curve_fit
from scipy.special import erfcx
from skimage import filters, io, morphology
from skimage.measure import label as skimage_label
from skimage.segmentation import expand_labels

from tttrkit.ptuio.decoder import T3OverflowCorrector
from tttrkit.ptuio.reconstructor import ScanConfig
from tttrkit.ptuio.utils import estimate_tcspc_bins

from flim_fret_data_processing.config import REPO_DIR

# ---------------------------------------------------------------------------
# Dataset layout
# ---------------------------------------------------------------------------

MASK_KINDS = ('ib', 'cell')


def index_samples(raw_dir):
    """{sample -> {'path', 'session', 'files'}} for `raw_dir/<session>/<sample>/*.ptu`.

    A sample is any session subfolder that holds .ptu files, so exported
    siblings (`<sample>_export`, `<sample>_bf`) are skipped. Sample names must
    be unique across sessions.
    """
    samples = {}
    for session in sorted(os.listdir(raw_dir)):
        session_dir = os.path.join(raw_dir, session)
        if not os.path.isdir(session_dir):
            continue
        for sample in sorted(os.listdir(session_dir)):
            sample_dir = os.path.join(session_dir, sample)
            if not os.path.isdir(sample_dir):
                continue
            files = sorted(f for f in os.listdir(sample_dir) if f.endswith('.ptu'))
            if not files:
                continue
            if sample in samples:
                raise ValueError(f"sample '{sample}' appears in both "
                                 f"{samples[sample]['session']} and {session}")
            samples[sample] = {'path': sample_dir, 'session': session, 'files': files}
    return samples


def index_masks(mask_root, samples):
    """{(sample, image base name, kind) -> mask path} for the curated masks.

    Mask folders are found by name at any depth under `mask_root`:
    `<sample>_ib_masks` and `<sample>_cell_masks` (an optional `_export` before
    the suffix is accepted). Inside, a mask belongs to the image whose base
    name it starts with. If an image has the same mask in several folders, the
    last folder in sorted order is used and the sample is reported.
    """
    mask_dirs = {}
    for dirpath, dirnames, _ in os.walk(mask_root):
        for d in dirnames:
            for kind in MASK_KINDS:
                suffix = f"_{kind}_masks"
                if not d.endswith(suffix):
                    continue
                sample = d[:-len(suffix)]
                if sample.endswith('_export'):
                    sample = sample[:-len('_export')]
                mask_dirs.setdefault((sample, kind), []).append(os.path.join(dirpath, d))

    masks, ambiguous = {}, set()
    for sample, info in samples.items():
        bases = [os.path.splitext(f)[0] for f in info['files']]
        for kind in MASK_KINDS:
            for mask_dir in sorted(mask_dirs.get((sample, kind), [])):
                for f in sorted(os.listdir(mask_dir)):
                    if not f.endswith('.tif'):
                        continue
                    for base in bases:
                        if f.startswith(base + '_'):
                            if (sample, base, kind) in masks:
                                ambiguous.add(sample)
                            masks[(sample, base, kind)] = os.path.join(mask_dir, f)

    if ambiguous:
        print(f"WARNING: several mask folders hold the same images for "
              f"{sorted(ambiguous)}; the last folder in sorted order is used.")
    return masks


def load_sample_rename(path):
    """{sample folder -> display name} from a two-column old_name/new_name CSV."""
    if not path:
        return {}
    key = pd.read_csv(path, dtype=str)
    key.columns = [c.strip() for c in key.columns]
    key = key[['old_name', 'new_name']].fillna('').apply(lambda col: col.str.strip())
    key = key[(key['old_name'] != '') & (key['new_name'] != '')]
    return dict(zip(key['old_name'], key['new_name']))


# ---------------------------------------------------------------------------
# Acquisition parameters and scan geometry
# ---------------------------------------------------------------------------

def read_acquisition_params(tags, target_bin_width_ps=40):
    """Acquisition parameters of one .ptu file, from its header tags.

    `tcspc_binning` is the factor that brings the raw TCSPC resolution to
    `target_bin_width_ps`, so that decays recorded at 5 ps and at 40 ps end up
    on the same time axis.
    """
    resolution = tags.get("MeasDesc_Resolution", 5e-12)
    n_raw_bins = estimate_tcspc_bins(tags, buffer=0)
    binning = max(1, round(target_bin_width_ps * 1e-12 / resolution))
    return {
        'wrap':             tags.get("TTResultFormat_WrapAround", 1024),
        'sync_rate':        tags.get("TTResult_SyncRate", 40e6),
        'tcspc_resolution': resolution,
        'tcspc_bins':       n_raw_bins,
        'tcspc_binning':    binning,
        'n_bins':           n_raw_bins // binning,
        'bin_width':        resolution * binning,
        'pixel_size_um':    tags.get("ImgHdr_PixResol"),
    }


def scan_shape(ptu_path, tags, frames=5, fallback_accumulations=4):
    """(lines, pixels, source) of one acquisition.

    The line length is in the .ptu header (`ImgHdr_PixX`), the image height is
    not: `ImgHdr_PixY` counts scanned lines, i.e. height x line accumulations,
    and no tag records the accumulations. The height is therefore read from the
    .msr file next to the .ptu; without one, `fallback_accumulations` is assumed.
    """
    pixels = tags['ImgHdr_PixX']
    scanned_lines = tags['ImgHdr_PixY']

    msr_path = os.path.splitext(ptu_path)[0] + '.msr'
    try:
        msr_pixels, lines = read_msr_scan_shape(msr_path, z_planes=frames)
        if msr_pixels != pixels:
            raise ValueError(f".msr says {msr_pixels} px per line, .ptu says {pixels}")
        if scanned_lines % lines:
            raise ValueError(f"{scanned_lines} scanned lines is not a whole "
                             f"number of {lines}-line frames")
        return lines, pixels, '.msr'
    except (OSError, KeyError, ValueError) as error:
        lines = scanned_lines // fallback_accumulations
        return lines, pixels, (f'no .msr, assumed {fallback_accumulations} '
                               f'accumulations ({error})')


def make_scan_config(tags, lines, pixels, frames=5):
    """ScanConfig of one acquisition on a (lines, pixels) grid."""
    scanned_lines = tags['ImgHdr_PixY']
    if scanned_lines % lines:
        raise ValueError(f"{scanned_lines} scanned lines is not a whole number "
                         f"of {lines}-line frames")
    return ScanConfig(
        bidirectional=False,
        frames=frames,
        lines=lines,
        pixels=pixels,
        line_accumulations=(scanned_lines // lines,),
        max_detector=2,
        line_start_marker_channel=(1, 3),   # line markers arrive in channels 1 and 3
    )


def run_reconstruction(reader, reconstructor, wraparound):
    """Feed all photon records of `reader` through `reconstructor`; returns its result.

    The reconstructor's console messages (one or more per file) are suppressed.
    """
    corrector = T3OverflowCorrector(wraparound=wraparound)
    with open(os.devnull, 'w') as devnull, contextlib.redirect_stdout(devnull):
        for chunk in reader.iter_chunks():
            reconstructor.update(corrector.correct(chunk))
        return reconstructor.finalize()


# ---------------------------------------------------------------------------
# Abberior .msr (Imspector OBF) reading
#
# The .ptu holds only the two photon-counting channels. The transmitted-light
# detector is analog, so the brightfield image exists only in the .msr written
# next to each .ptu. It is recorded on the same scan as the FLIM z-stack and
# normally shares its pixel grid; `check_bf_registration` catches the files
# where it does not. The OBF container is parsed directly: a chain of stacks,
# each a fixed-size header followed by a zlib-compressed data blob.
# ---------------------------------------------------------------------------

OBF_STACK_MAGIC = b'OMAS_BF_STACK\n\xff\xff'

# OBF value types are bit flags; this microscope writes 0x08 (int16).
OBF_DTYPES = {
    0x01: 'u1', 0x02: 'i1',
    0x04: 'u2', 0x08: 'i2',
    0x10: 'u4', 0x20: 'i4',
    0x40: 'f4', 0x80: 'f8',
}


def read_msr_stacks(msr_path, load=True):
    """List every stack in an Abberior .msr file.

    Returns one dict per stack with:
        'name'  channel name as shown in Imspector, e.g. 'BF {9}'
        'shape' (x, y, z, t) sizes, in Imspector's axis order
        'phys'  physical extent of each axis, in metres
        'off'   physical offset of each axis, in metres (stage position)
        'data'  numpy array in (t, z, y, x) order -- only if load=True
    """
    with open(msr_path, 'rb') as fh:
        raw_file = fh.read()

    stacks = []
    search_from = 0

    while True:
        magic_at = raw_file.find(OBF_STACK_MAGIC, search_from)
        if magic_at < 0:
            break
        search_from = magic_at + 1

        p = magic_at + len(OBF_STACK_MAGIC)
        version, rank = struct.unpack_from('<2i', raw_file, p); p += 8
        sizes = struct.unpack_from('<15i', raw_file, p); p += 60
        lengths = struct.unpack_from('<15d', raw_file, p); p += 120
        offsets = struct.unpack_from('<15d', raw_file, p); p += 120
        dtype_flag, compression, _level = struct.unpack_from('<3i', raw_file, p); p += 12
        name_len, desc_len = struct.unpack_from('<2i', raw_file, p); p += 8
        p += 8  # reserved
        data_len, _next_stack = struct.unpack_from('<2q', raw_file, p); p += 16
        name = raw_file[p:p + name_len].decode('utf-8', 'replace'); p += name_len
        description = raw_file[p:p + desc_len].decode('utf-8', 'replace'); p += desc_len

        shape = tuple(sizes[:rank])
        stack = {
            'name': name,
            'description': description,
            'shape': shape,
            'phys': tuple(lengths[:rank]),
            'off': tuple(offsets[:rank]),
            'version': version,
            'offset': magic_at,
        }

        if load:
            blob = raw_file[p:p + data_len]
            buf = zlib.decompress(blob) if compression else blob

            if dtype_flag not in OBF_DTYPES:
                raise ValueError(
                    f"{os.path.basename(msr_path)}: stack {name!r} has unknown OBF "
                    f"value type {dtype_flag}"
                )
            dtype = np.dtype('<' + OBF_DTYPES[dtype_flag])

            expected = int(np.prod(shape)) * dtype.itemsize
            if len(buf) != expected:
                raise ValueError(
                    f"{os.path.basename(msr_path)}: stack {name!r} decoded to "
                    f"{len(buf)} bytes, expected {expected} for shape {shape} "
                    f"of {dtype}"
                )

            # Imspector stores x fastest, so the buffer is (t, z, y, x).
            stack['data'] = np.frombuffer(buf, dtype=dtype).reshape(shape[::-1])

        stacks.append(stack)

    return stacks


def find_msr_flim_stack(msr_path, z_planes=5, exclude_names=('BF', 'Ch3')):
    """Return the .msr stack describing the FLIM acquisition, as a stack dict.

    One .msr holds an overview, preview snapshots, the transmission channel and
    the fluorescence z-stack the .ptu photons belong to. The last is the only
    fluorescence stack with `z_planes` planes; of the two detector channels
    that qualify, the larger one is returned.

    Raises KeyError if no fluorescence stack has that depth.
    """
    candidates = []
    for stack in read_msr_stacks(msr_path, load=False):
        name = stack['name'].split(' {')[0]
        if name in exclude_names or len(stack['shape']) < 3:
            continue
        if stack['shape'][2] == z_planes:  # 'shape' is Imspector's (x, y, z, t)
            candidates.append(stack)

    if not candidates:
        raise KeyError(
            f"{os.path.basename(msr_path)}: no fluorescence stack with "
            f"{z_planes} z planes"
        )

    return max(candidates, key=lambda s: s['shape'][0] * s['shape'][1])


def read_msr_scan_shape(msr_path, z_planes=5):
    """(pixels, lines) of the FLIM scan, i.e. its true width and height."""
    shape = find_msr_flim_stack(msr_path, z_planes)['shape']
    return int(shape[0]), int(shape[1])


def check_bf_registration(msr_path, z_planes=5, tolerance_px=0.5):
    """Check that the transmission stack covers the FLIM field of view.

    Returns (ok, reason). Imspector sometimes leaves the brightfield stack of
    the previous acquisition in the file: a sharp image of a different field
    of view, usually with the same shape. The stage offsets of the two stacks
    are therefore compared, in units of FLIM pixels.
    """
    try:
        flim = find_msr_flim_stack(msr_path, z_planes)
    except KeyError as e:
        return False, str(e)

    bf_stacks = [s for s in read_msr_stacks(msr_path, load=False)
                 if s['name'].split(' {')[0] in ('BF', 'Ch3')]
    if not bf_stacks:
        return False, 'no transmission channel'
    bf = max(bf_stacks, key=lambda s: s['shape'][2] if len(s['shape']) > 2 else 1)

    # Pixel size of the FLIM scan, per axis, in metres.
    pixel_size = [flim['phys'][a] / flim['shape'][a] for a in (0, 1)]

    shift_px = [(bf['off'][a] - flim['off'][a]) / pixel_size[a] for a in (0, 1)]
    span_px = [(bf['phys'][a] - flim['phys'][a]) / pixel_size[a] for a in (0, 1)]

    def fmt(values):
        return '(' + ', '.join(f'{v:+.1f}' for v in values) + ')'

    if max(abs(s) for s in shift_px) > tolerance_px:
        magnitude = ('stale BF stack -- it holds a different field of view'
                     if max(abs(s) for s in shift_px) > 5
                     else 'BF and FLIM drifted apart')
        return False, (f"BF is offset from FLIM by {fmt(shift_px)} px in xy "
                       f"-- {magnitude}")

    if max(abs(s) for s in span_px) > tolerance_px:
        return False, (f"BF spans {fmt(span_px)} px more than FLIM in xy "
                       f"-- different field size")

    if bf['shape'][:2] != flim['shape'][:2]:
        return False, (f"BF is {bf['shape'][:2]} px, FLIM is "
                       f"{flim['shape'][:2]} px -- different sampling")

    return True, 'BF and FLIM share the same grid'


def read_msr_transmission(msr_path, channel_names=('BF', 'Ch3')):
    """Return the transmitted-light z-stack of an .msr as a (z, y, x) array.

    `channel_names` are tried in order, with the Imspector config suffix
    stripped ('BF {9}' -> 'BF'). Raises KeyError if the file has none of them.
    """
    stacks = read_msr_stacks(msr_path, load=True)
    by_name = {}
    for stack in stacks:
        by_name.setdefault(stack['name'].split(' {')[0], []).append(stack)

    for wanted in channel_names:
        if wanted not in by_name:
            continue
        # Prefer the deepest stack, so a z-stack beats a single-plane snapshot.
        stack = max(by_name[wanted],
                    key=lambda s: s['shape'][2] if len(s['shape']) > 2 else 1)
        data = stack['data']
        # (t, z, y, x) -> (z, y, x); these are single-timepoint acquisitions.
        if data.ndim == 4:
            if data.shape[0] != 1:
                raise ValueError(
                    f"{os.path.basename(msr_path)}: {stack['name']!r} has "
                    f"{data.shape[0]} timepoints, expected 1"
                )
            data = data[0]
        elif data.ndim == 2:  # a single-plane channel, e.g. 'Ch3'
            data = data[np.newaxis]
        return np.ascontiguousarray(data)

    raise KeyError(
        f"{os.path.basename(msr_path)}: no transmission channel "
        f"(tried {', '.join(channel_names)}); found "
        f"{', '.join(sorted(by_name)) or 'no stacks at all'}"
    )


# ---------------------------------------------------------------------------
# Segmentation (notebooks 02a-c)
# ---------------------------------------------------------------------------

def segment_ibs(image, sigma_small=1, sigma_large=3, min_size=10, expand_distance=1.8):
    """One label per inclusion body, from a difference-of-Gaussians band-pass.

    The DoG is computed in raw photon counts and thresholded by Otsu's method
    over its positive values only, i.e. over the pixels brighter than their
    surroundings. Objects below `min_size` px are dropped and every label is
    grown by `expand_distance` px to take in the IB's edge pixels.
    """
    dog = (filters.gaussian(image, sigma=sigma_small, preserve_range=True)
           - filters.gaussian(image, sigma=sigma_large, preserve_range=True))

    ib_mask = dog > filters.threshold_otsu(dog[dog > 0])
    ib_mask = morphology.remove_small_objects(ib_mask, min_size=min_size)

    return expand_labels(skimage_label(ib_mask), distance=expand_distance)


def keep_solid_pieces(labels, min_area, min_half_thickness):
    """Split every label into connected pieces and keep only the solid ones.

    Cutting the IB halo out of a cell mask leaves stray pixels and thin rinds
    along the halo boundary, which carry IB bleed-through rather than cytosol.
    Each connected piece is judged on its own area and on its half thickness
    (the radius of the largest disk that fits inside it); survivors are
    relabelled 1..n.

    Returns the surviving pieces and the number dropped.
    """
    kept = np.zeros_like(labels)
    n_kept = n_dropped = 0

    for lbl in range(1, int(labels.max()) + 1):
        pieces, n_pieces = ndimage.label(labels == lbl)
        for piece_id in range(1, n_pieces + 1):
            piece = pieces == piece_id
            half_thickness = ndimage.distance_transform_edt(piece).max()
            if piece.sum() < min_area or half_thickness < min_half_thickness:
                n_dropped += 1
                continue
            n_kept += 1
            kept[piece] = n_kept

    return kept, n_dropped


def normalize_for_display(image, low_percentile=1, high_percentile=99.5):
    """Percentile-stretch an image to 0-1 for display."""
    lo, hi = np.percentile(image, [low_percentile, high_percentile])
    if hi <= lo:
        return np.zeros_like(image, dtype=np.float32)
    return np.clip((image.astype(np.float32) - lo) / (hi - lo), 0, 1)


def overlay_labels(image, labels, alpha=0.6, seed=0):
    """Alpha-blend coloured label regions onto a percentile-stretched image."""
    rgb = np.stack([normalize_for_display(image)] * 3, axis=-1)

    n_labels = int(labels.max())
    if n_labels == 0:
        return (rgb * 255).astype(np.uint8)

    cmap = plt.get_cmap('gist_rainbow', n_labels)
    colors = np.array([cmap(i)[:3] for i in range(n_labels)])
    colors = colors[np.random.default_rng(seed).permutation(n_labels)]

    for i in range(n_labels):
        region = labels == (i + 1)
        rgb[region] = (1 - alpha) * rgb[region] + alpha * colors[i]

    return (rgb * 255).astype(np.uint8)


def save_labels(folder, file_name, labels):
    """Write a label image as uint8 to `folder/file_name`."""
    if labels.max() > np.iinfo(np.uint8).max:
        raise ValueError(f"{file_name}: {labels.max()} labels do not fit in uint8")
    os.makedirs(folder, exist_ok=True)
    io.imsave(os.path.join(folder, file_name), labels.astype(np.uint8),
              check_contrast=False)


def save_overview(path, title, rows, column_titles):
    """One .png per sample: a row of overlays per image.

    `rows` is a list of (image name, overlay, overlay, ...).
    """
    save_path = Path(path)
    n_columns = len(column_titles)
    fig, axes = plt.subplots(len(rows), n_columns,
                             figsize=(3 * n_columns, 3 * len(rows)), squeeze=False)

    for row, (file_name, *overlays) in enumerate(rows):
        for col, overlay in enumerate(overlays):
            axes[row, col].imshow(overlay)
            axes[row, col].set_xticks([])
            axes[row, col].set_yticks([])
        axes[row, 0].set_ylabel(file_name, fontsize=8)

    for col, column_title in enumerate(column_titles):
        axes[0, col].set_title(column_title)

    fig.suptitle(title)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(save_path, dpi=150)
    plt.close(fig)
    print(f"Overview saved to {save_path.relative_to(REPO_DIR)}")


def show_image_grid(images, title, cols=4, tile_size=(2.5, 2.5)):
    """Show a list of images in a grid."""
    rows = (len(images) + cols - 1) // cols
    fig, axes = plt.subplots(rows, cols, squeeze=False,
                             figsize=(cols * tile_size[0], rows * tile_size[1]))
    for ax in axes.ravel():
        ax.axis('off')
    for ax, image in zip(axes.ravel(), images):
        ax.imshow(image)
    fig.suptitle(title)
    plt.tight_layout()
    plt.show()


# ---------------------------------------------------------------------------
# Phasors
#
# The phasor used throughout is <exp(+i*omega*t)>, so a single-exponential
# decay sits at 1 / (1 - i*omega*tau), with s > 0.
# ---------------------------------------------------------------------------

def calculate_phasor_from_lifetime(lifetime, frequency):
    """Phasor of a single-exponential decay (lifetime in s, frequency in Hz)."""
    return 1 / (1 - 1j * 2 * np.pi * frequency * lifetime)


def apparent_lifetimes(phasor, sync_rate):
    """(tau_phi, tau_m) in seconds: phase and modulation lifetime of a phasor."""
    z = np.asarray(phasor, dtype=np.complex128)
    omega = 2.0 * np.pi * sync_rate
    with np.errstate(divide='ignore', invalid='ignore'):
        tau_phi = np.imag(z) / (omega * np.real(z))
        modulation_sq = np.abs(z) ** 2
        tau_m = np.sqrt(np.clip(1.0 / modulation_sq - 1.0, 0.0, None)) / omega
    return tau_phi, tau_m


# ---------------------------------------------------------------------------
# Leading-edge timing
#
# Between imaging sessions the excitation pulse shifts within the TCSPC window.
# A histogram delayed by dt has its phasor multiplied by exp(+i*omega*dt): a
# pure rotation, undone by exp(-i*omega*dt). The delay is read from the rising
# edge, which below ~30% of the peak is dominated by the instrument response
# and is nearly independent of the fluorescence lifetime.
# ---------------------------------------------------------------------------

def rise_baseline(decay, smoothing=1.0):
    """(baseline, peak index, peak height above baseline) of a TCSPC histogram.

    The baseline is the median of the bins before the rise starts; the start of
    the rise is located from a first percentile-based guess.
    """
    y = gaussian_filter1d(np.asarray(decay, dtype=float), smoothing) if smoothing else \
        np.asarray(decay, dtype=float)

    peak_index = int(np.argmax(y))
    if peak_index < 4:
        raise ValueError("peak sits at the very start of the histogram -- no rising edge to fit")

    guess = np.percentile(y[:peak_index], 10)
    above = np.flatnonzero(y[:peak_index] - guess < 0.01 * (y[peak_index] - guess))
    rise_start = int(above[-1]) if len(above) else peak_index // 2

    baseline = float(np.median(y[:max(rise_start, 3)]))
    return baseline, peak_index, float(y[peak_index] - baseline)


def leading_edge_crossing(decay, bin_width, fractions=(0.05, 0.10, 0.15, 0.20, 0.25),
                          smoothing=1.0, return_spread=False):
    """Timing offset of a TCSPC histogram, from constant-fraction crossings of its rise.

    The rise is baseline-subtracted and normalised to the peak; the times at
    which it crosses each of `fractions` are linearly interpolated and averaged.

    Parameters
    ----------
    decay : array_like
        Accumulated TCSPC histogram.
    bin_width : float
        Width of one channel, in seconds.
    fractions : sequence of float
        Fractions of the peak height to cross. Keep them below ~0.3.
    smoothing : float
        Gaussian sigma (channels) used only to find the baseline and the peak.
    return_spread : bool
        Also return the standard deviation across the individual crossings.

    Returns
    -------
    float, or (float, float) when `return_spread`
        Offset in seconds. Only differences between measurements are meaningful.
    """
    baseline, peak_index, height = rise_baseline(decay, smoothing)
    if height <= 0:
        return (np.nan, np.nan) if return_spread else np.nan

    y = (np.asarray(decay, dtype=float)[:peak_index + 1] - baseline) / height

    crossings = []
    for fraction in np.atleast_1d(fractions):
        below = np.flatnonzero(y < fraction)
        if not len(below) or below[-1] + 1 > peak_index:
            continue
        i = int(below[-1])
        step = y[i + 1] - y[i]
        if step <= 0:
            continue
        crossings.append(i + (fraction - y[i]) / step)

    if not crossings:
        return (np.nan, np.nan) if return_spread else np.nan

    offset = float(np.mean(crossings)) * bin_width
    if return_spread:
        return offset, float(np.std(crossings)) * bin_width
    return offset


def _emg_model(t, t0, sigma, tau, amplitude, baseline):
    """Gaussian IRF convolved with a single exponential (exponentially modified Gaussian)."""
    u = t - t0
    z = (sigma ** 2 / tau - u) / (np.sqrt(2.0) * sigma)
    return baseline + amplitude * 0.5 * erfcx(z) * np.exp(-0.5 * (u / sigma) ** 2)


def fit_leading_edge(decay, bin_width, lo=0.03, hi=0.60, tail_bins=25, smoothing=1.0):
    """Timing offset and IRF width of a TCSPC histogram, fitted to its rising edge.

    Model-based cross-check of `leading_edge_crossing`: an exponentially
    modified Gaussian is fitted from the `lo` crossing of the rise to
    `tail_bins` channels past the `hi` crossing, weighted for Poisson noise.

    Returns a dict with t0 and sigma (Gaussian IRF width) in seconds, tau (the
    nuisance lifetime), t0_err, the photons in the fit window and the window.
    """
    counts = np.asarray(decay, dtype=float)
    baseline, peak_index, height = rise_baseline(counts, smoothing)
    if height <= 0:
        raise ValueError("histogram has no peak above its baseline")

    smooth = gaussian_filter1d(counts, smoothing) if smoothing else counts
    normalised = (smooth[:peak_index + 1] - baseline) / height

    below_lo = np.flatnonzero(normalised < lo)
    below_hi = np.flatnonzero(normalised < hi)
    if not len(below_lo) or not len(below_hi):
        raise ValueError("could not locate the rise window")

    first = int(below_lo[-1])
    last = min(int(below_hi[-1]) + 1 + tail_bins, len(counts) - 1)
    if last - first < 6:
        raise ValueError(f"rise window is only {last - first} channels wide")

    index = np.arange(first, last + 1, dtype=float)
    window_counts = counts[first:last + 1]
    weights = np.sqrt(window_counts + 1.0)

    # Fitted in channels, so all parameters are of similar magnitude.
    sigma_guess = max(1.0, 0.5 * (last - first))
    t0_guess = 0.5 * (first + last)
    p0 = [t0_guess, sigma_guess, 10 * sigma_guess, height, baseline]
    bounds = ([first - 5 * sigma_guess, 1e-3, 1e-2, 0.0, -abs(height)],
              [last + 5 * sigma_guess, 10 * sigma_guess, 1e4, 1e3 * height, height])

    popt, pcov = curve_fit(_emg_model, index, window_counts, p0=p0, sigma=weights,
                           absolute_sigma=False, bounds=bounds, maxfev=20000)

    return {
        't0':      popt[0] * bin_width,
        't0_err':  float(np.sqrt(pcov[0, 0])) * bin_width,
        'sigma':   popt[1] * bin_width,
        'tau':     popt[2] * bin_width,
        'photons': float(window_counts.sum()),
        'window':  (first, last),
    }


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

# Constructs are named <family>-<number> (H2-5, L1-3); the family is the
# selection round and has one colour in every figure.
FAMILY_COLORS = {
    'H1': '#79d2d7ff',
    'H2': '#2f8cbfff',
    'H3': '#1446a1ff',
    'L1': '#a66a20ff',
    'L2': '#7e1700ff',
}

# Reference constructs and controls, which belong to no family.
FAMILY_COLOR_FALLBACK = '#6e6e6eff'


def family_of(display_name):
    """Family of a construct display name ('H2-5' -> 'H2'), else 'other'."""
    match = re.match(r'^([A-Za-z]+\d*)', str(display_name))
    return match.group(1) if match and match.group(1) in FAMILY_COLORS else 'other'


def family_color(name):
    """Colour of a family ('H2') or of a construct display name ('H2-5')."""
    return FAMILY_COLORS.get(name) or FAMILY_COLORS.get(family_of(name), FAMILY_COLOR_FALLBACK)


def save_figure(fig, name, directory, formats=('svg',), dpi=300, **kwargs):
    """Write `fig` to `directory/name.<ext>` for each format; returns the paths."""
    os.makedirs(directory, exist_ok=True)

    stem = "".join(character if character.isalnum() or character in "-_." else "_"
                   for character in str(name))

    paths = []
    for fmt in formats:
        path = os.path.join(directory, f"{stem}.{fmt}")
        raster = {} if fmt in ('svg', 'pdf', 'eps', 'ps') else {'dpi': dpi}
        fig.savefig(path, format=fmt, bbox_inches='tight', **raster, **kwargs)
        paths.append(path)

    return paths


def draw_universal_circle(ax, sync_rate, tau_max=None, tick_length=0.02,
                          color='black', label_color='black', xlim=None, ylim=None):
    """Draw the universal semicircle with lifetime ticks (in ns) on `ax`.

    Ticks and labels that fall outside `xlim`/`ylim` are left out, so the same
    call serves the full plane and a zoomed window.
    """
    if xlim is None:
        xlim = (-0.1, 1.1)
    if ylim is None:
        ylim = (0, 0.8)

    def inside_limits(point):
        x, y = point
        return xlim[0] <= x <= xlim[1] and ylim[0] <= y <= ylim[1]

    omega = 2 * np.pi * sync_rate

    if tau_max is None:
        tau_max = int(np.ceil(1 / sync_rate / 2))

    ticks = np.arange(1, tau_max + 1)
    taus_ns = ticks[(ticks <= 8) | ((ticks > 8) & (ticks % 2 == 0))]

    center = np.array([0.5, 0])
    radius = 0.5

    theta = np.linspace(0, np.pi, 300)
    ax.plot(center[0] + radius * np.cos(theta), center[1] + radius * np.sin(theta),
            '-', color=color, label='Universal Circle', lw=1)

    for tau_ns in taus_ns:
        tau = tau_ns * 1e-9
        p = np.array([1 / (1 + (omega * tau) ** 2),
                      (omega * tau) / (1 + (omega * tau) ** 2)])

        v = p - center
        v_unit = v / np.linalg.norm(v)

        p1 = p - (tick_length / 2) * v_unit
        p2 = p + (tick_length / 2) * v_unit
        label_pos = p + (tick_length * 1.2) * v_unit

        if inside_limits(p1) or inside_limits(p2):
            ax.plot([p1[0], p2[0]], [p1[1], p2[1]], '-', color=color, lw=1)

        if inside_limits(label_pos):
            ax.text(label_pos[0], label_pos[1], f'{tau_ns}', fontsize=8,
                    ha='center', va='center', color=label_color)

    ax.set_xlabel('g')
    ax.set_ylabel('s')
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_aspect('equal', adjustable='datalim')


# sRGB (D65) <-> CIE Lab, for the isoluminant colourmap of the lifetime maps.
_XYZ_FROM_RGB = np.array([[0.4124564, 0.3575761, 0.1804375],
                          [0.2126729, 0.7151522, 0.0721750],
                          [0.0193339, 0.1191920, 0.9503041]])
_RGB_FROM_XYZ = np.linalg.inv(_XYZ_FROM_RGB)
_WHITE_POINT = np.array([0.95047, 1.0, 1.08883])
_DELTA = 6 / 29


def _lab_to_linear_rgb(lightness, a, b):
    """Linear-light sRGB of a Lab colour. Values outside [0, 1] are out of gamut."""
    lightness, a, b = np.broadcast_arrays(*(np.asarray(v, float) for v in (lightness, a, b)))
    fy = (lightness + 16) / 116
    f = np.stack([fy + a / 500, fy, fy - b / 200], axis=-1)
    xyz = np.where(f > _DELTA, f ** 3, 3 * _DELTA ** 2 * (f - 4 / 29)) * _WHITE_POINT
    return xyz @ _RGB_FROM_XYZ.T


def _linear_rgb_to_srgb(linear):
    linear = np.clip(linear, 0, 1)
    return np.where(linear <= 0.0031308, 12.92 * linear,
                    1.055 * linear ** (1 / 2.4) - 0.055)


def isoluminant_colormap(lightness=65.0, hue_range=(270, 30), n=256, name='isoluminant'):
    """Hue sweep at constant CIE L*, at the largest chroma inside the sRGB gamut.

    Constant lightness keeps the hue independent of brightness, so a second
    quantity can be encoded in the brightness of the same pixel. `hue_range`
    is in CIELCh degrees. Returns (colormap, chroma).
    """
    hue = np.radians(np.linspace(*hue_range, n))
    chroma = 0.0
    for candidate in np.arange(70, 2, -0.5):
        linear = _lab_to_linear_rgb(lightness,
                                    candidate * np.cos(hue), candidate * np.sin(hue))
        if linear.min() >= -1e-4 and linear.max() <= 1 + 1e-4:
            chroma = candidate
            break
    if not chroma:
        raise RuntimeError(f"no in-gamut chroma at L*={lightness} over hues {hue_range}")

    colors = _linear_rgb_to_srgb(
        _lab_to_linear_rgb(lightness, chroma * np.cos(hue), chroma * np.sin(hue)))
    return matplotlib.colors.ListedColormap(colors, name=name), chroma
