# FLIM phasor analysis of inclusion bodies

Image processing and FLIM analysis code for **[paper title, authors, journal, year — add the citation]**.

The notebooks reconstruct intensity images from time-tagged photon data
(PicoQuant `.ptu`), segment inclusion bodies (IBs) and cytosol, and quantify
the donor (mTurquoise2) lifetime and apparent FRET efficiency in both
compartments by phasor analysis.

The raw data, curated masks and result tables are available at
**[Zenodo DOI — add the link]**.

## Workflow

| Step | File | What it does |
|---|---|---|
| 1 | `01_export_intensity.ipynb` | Intensity images from the `.ptu` files and brightfield stacks from the `.msr` files |
| 2 | `02a_segmentation_brightfield.ipynb` | IBs from fluorescence, cells from brightfield (Cellpose-SAM) |
| 2 | `02b_segmentation_fluorescence.ipynb` | IBs and cells from fluorescence, for images without usable brightfield |
| 2 | `02c_segmentation_diffuse.ipynb` | Whole cells, for samples without IBs |
| 3 | `curate_masks.py` | Manual inspection and correction of the masks in napari |
| 4 | `03_phasor_analysis.ipynb` | Per-mask phasors, session correction, donor calibration; tables of τ<sub>φ</sub>, τ<sub>m</sub>, E<sub>φ</sub>, E<sub>m</sub> |
| 5 | `04_lifetime_maps.ipynb` | Intensity-modulated lifetime maps, per image and per construct |
| 5 | `05_construct_plots.ipynb` | Phasor and FRET efficiency plots per construct |

`helpers.py` holds the functions the notebooks share and `config.py` the data
locations. `rename_key.csv` maps sample folder names to the construct names
used in the tables and figures.

Each sample goes through one of the three segmentation notebooks; the sample
lists at the top of each notebook record which. Segmentation parameters were
tuned between batches, and the notebooks hold the values of the final runs.
All masks were curated by hand afterwards, so the curated masks in the dataset,
not a fresh segmentation, are the reference for the analysis.

Notebooks 03-05 only need the raw data and the curated masks, so the published
results can be reproduced without repeating steps 1-3. The whole dataset
(271 images) takes a few minutes per notebook on a laptop CPU.

## Installation

```bash
conda env create -f requirements.yml
conda activate flim-buchel-2025
jupyter lab
```

Photon data are read with [tttrkit](https://github.com/IMCF-Biocev/tttrkit).
Tested with Python 3.10, tttrkit 0.1.0, scikit-image 0.25, Cellpose 4.2 and
napari 0.7 on Linux; a GPU is not required.

## Data layout

Unpack the dataset into `data/` next to the notebooks, or set `DATA_ROOT` in
`config.py` (or the `FLIM_DATA_ROOT` environment variable) to its location.

```
data/
├── raw/                          acquisitions, one folder per imaging session
│   └── 2025-07-30/
│       └── 12/                   one folder per sample (construct)
│           ├── <image>.ptu       time-tagged photon data
│           └── <image>.msr       Imspector file with the brightfield stack (if recorded)
├── intensity/                    written by notebook 01
│   └── 2025-07-30/
│       ├── 12_export/<image>.tif       intensity image
│       └── 12_bf/<image>.tif           brightfield z-stack
├── masks_auto/                   written by notebooks 02a-c
│   ├── 12_ib_masks/<image>_ib_labels.tif
│   ├── 12_cell_masks/<image>_cytosol_labels.tif
│   └── 12_bf_overview.png              overlays for a quick check
├── masks_curated/                masks after manual curation, read by notebooks 03 and 04
│   ├── 12_ib_masks/...
│   └── 12_cell_masks/...
└── results/                      written by notebooks 03-05
    ├── phasors.csv, decays.csv, fret_objects.csv, fret_images.csv, fret_constructs.csv
    ├── figures/
    └── lifetime_maps/
```

Sample names must be unique across sessions. Mask folders are recognised by
name at any depth below `masks_curated/` (`<sample>_ib_masks`,
`<sample>_cell_masks`, optionally `<sample>_export_...`), and a mask belongs to
the image whose name it starts with. Samples without IBs have only a cell
mask (`<image>_cell_labels.tif`).

`masks_auto/` in the dataset holds the automatic masks as they were before
curation. It covers more images than `masks_curated/`, because images that were
discarded during curation have automatic masks only. Notebooks 02a-c write to
`masks_auto/`; to re-run the segmentation without replacing the deposited
masks, set `output_dir` in those notebooks to another folder.

## Mask curation

```bash
python curate_masks.py
```

Select the folder with the intensity images of a sample (`<sample>_export`),
its label folders (`<sample>_ib_masks`, `<sample>_cell_masks`; Cancel when
done) and the output folder (`masks_curated`). Each image opens in napari with
its label layers. Edit them with the label tools, then press `s` to save and
go to the next image, `n` to skip, `p` to go back. To discard an image, delete
all its labels and save.

## Results

| File | Content |
|---|---|
| `phasors.csv` | One row per mask file: raw and calibrated phasor, τ<sub>φ</sub>, τ<sub>m</sub> |
| `decays.csv` | Accumulated decay per sample and compartment, binned to 40 ps |
| `fret_objects.csv` | One row per labelled object (IB or cytosolic region) |
| `fret_images.csv` | Objects pooled per image |
| `fret_constructs.csv` | Objects pooled per construct, with the spread over objects and images |
| `lifetime_maps/` | Per-image panels and per-construct summaries of the lifetime maps |

Pooled values are photon-weighted: the phasor of a group is the sum of the
per-pixel phasor sums over the sum of the photon counts, which equals the
phasor of the summed decay.
