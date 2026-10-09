# FLIM phasor analysis of inclusion bodies

Image processing and FLIM analysis code for **[paper title, authors, journal, year — add the citation]**.

The notebooks reconstruct intensity images from time-tagged photon data
(PicoQuant `.ptu`), segment inclusion bodies (IBs) and cytosol, and quantify
the donor (mTurquoise2) lifetime and apparent FRET efficiency in both
compartments by phasor analysis.

The raw data, curated masks and result tables are available at
**[ZENODO](https://zenodo.org/records/23242459)**.

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
git clone https://github.com/buchelf-cuni/flim-fret-data-processing
cd flim-fret-data-processing
uv sync
```

Raw photon data are read with [tttrkit](https://github.com/IMCF-Biocev/tttrkit).
Tested with Python 3.10, tttrkit 0.1.0, scikit-image 0.25, Cellpose 4.2 and
napari 0.7 on Linux; a GPU is not required.

## Data layout

See `data\README.md` for details regarding the location of the data for this project.

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
