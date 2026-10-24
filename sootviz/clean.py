"""Cleaning a soot-foil record before anything is measured on it.

Three steps, in this order, each at a scale tied to the expected cell width:

  despeckle   grey opening then closing, to remove dust specks and pits that
              would otherwise be detected as structure;
  equalise    subtract the local mean and divide by the local standard
              deviation over a window SEVERAL CELLS WIDE. This is the step
              that matters on photographed foils, where illumination varies
              across the plate. The window must be well above the cell width:
              a window comparable to a cell removes the pattern itself;
  denoise     a light Gaussian well below the cell width.
"""
from __future__ import annotations

import numpy as np
from scipy import ndimage as ndi

from .filters import blur


def despeckle(img: np.ndarray, radius_px: int) -> np.ndarray:
    if radius_px < 1:
        return img
    r = int(radius_px)
    yy, xx = np.ogrid[-r:r + 1, -r:r + 1]
    disc = (xx ** 2 + yy ** 2) <= r * r
    return ndi.grey_closing(ndi.grey_opening(img, footprint=disc), footprint=disc)


def equalise(img: np.ndarray, window_px: float, clip: float = 3.0) -> np.ndarray:
    """Local mean/standard-deviation normalisation (a CLAHE stand-in)."""
    mu = blur(img, window_px)
    var = blur(img * img, window_px) - mu * mu
    sd = np.sqrt(np.clip(var, 0, None)) + 1e-6
    return np.clip((img - mu) / sd, -clip, clip)


def clean(img: np.ndarray, s: float, cell: float,
          despeckle_frac: float = 0.03,
          equalise_frac: float = 2.5,
          denoise_frac: float = 0.05) -> np.ndarray:
    """Full cleanup. Scales are fractions of the expected cell width `cell`."""
    out = despeckle(img, int(round(despeckle_frac * cell * s)))
    out = equalise(out, equalise_frac * cell * s)
    return ndi.gaussian_filter(out, denoise_frac * cell * s)
