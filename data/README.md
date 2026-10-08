# Where to find data

Download from ZENODO and copy to this folder

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
