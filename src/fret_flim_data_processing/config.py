"""Data locations shared by all notebooks.

Unpack the dataset into `data/` next to this file, or point `DATA_ROOT` (or
the FLIM_DATA_ROOT environment variable) at it. See the README for the layout.
"""
import os
from pathlib import Path

REPO_DIR = Path(__file__).resolve().parent

DATA_ROOT = Path(os.environ.get('FLIM_DATA_ROOT', REPO_DIR / 'data'))

# Raw acquisitions: <session>/<sample>/*.ptu (+ the matching Imspector .msr).
RAW_DIR = DATA_ROOT / 'raw'

# Intensity and brightfield images written by notebook 01:
# <session>/<sample>_export/*.tif and <session>/<sample>_bf/*.tif.
INTENSITY_DIR = DATA_ROOT / 'intensity'

# Masks written by notebooks 02a-c: <sample>_ib_masks/ and <sample>_cell_masks/.
MASKS_AUTO_DIR = DATA_ROOT / 'masks_auto'

# Masks after manual curation (curate_masks.py); read by notebooks 03 and 04.
MASKS_DIR = DATA_ROOT / 'masks_curated'

# Tables and figures written by notebooks 03-05.
RESULTS_DIR = Path(os.environ.get('FLIM_RESULTS_DIR', DATA_ROOT / 'results'))

# Sample folder name -> construct name used in tables and figures.
SAMPLE_NAMES_CSV = REPO_DIR / 'rename_key.csv'
