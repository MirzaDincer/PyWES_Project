"""Central project paths, so notebooks and scripts don't hard-code folders.

Usage:
    from pywes.paths import RAW_DIR, PROCESSED_DIR, FIGURES_DIR
"""

from pathlib import Path

# Repository root = two levels above this file (src/pywes/paths.py)
ROOT_DIR = Path(__file__).resolve().parents[2]

DATA_DIR = ROOT_DIR / "data"
RAW_DIR = DATA_DIR / "raw"
PROCESSED_DIR = DATA_DIR / "processed"
SAMPLE_DIR = DATA_DIR / "sample"

RESULTS_DIR = ROOT_DIR / "results"
FIGURES_DIR = RESULTS_DIR / "figures"
