"""Thin wrapper around the original observation pipeline."""
from _process_data import get_arrays
from .config import FIT_START, FIT_END
def fit_arrays():
    """Return the original 343-observation fit window."""
    return get_arrays(FIT_START, FIT_END)
