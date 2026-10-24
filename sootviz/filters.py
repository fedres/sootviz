"""Smoothing whose cost does not depend on the width of the kernel.

Every length in this pipeline is a fraction of the cell width, and the cell
width is what is being searched for, so the search evaluates the pipeline at
coarse trial scales where a Gaussian would need a kernel hundreds of pixels
wide. Three passes of a box filter have the same second moment as a Gaussian
and cost the same per pixel whatever that moment is, because a box filter is a
running sum. For a box of width w the variance is (w^2 - 1)/12, so three
passes match a Gaussian of width sigma at w = sqrt(4 sigma^2 + 1).

The approximation is used only where the kernel is wide, which is also where
its exact shape does not matter: the along-track smoothing of a matched filter
and the window of the local contrast equalisation. Narrow kernels -- the
across-track width of the matched filter, the fine scale of the ridge map --
stay exact, because there the shape is the whole point.
"""
from __future__ import annotations

import numpy as np
from scipy import ndimage as ndi

EXACT_BELOW = 2.0       # sigma in px below which the exact Gaussian is cheap


def blur1d(a: np.ndarray, sigma: float, axis: int) -> np.ndarray:
    """Smooth along one axis by sigma pixels, in O(1) per pixel."""
    if sigma <= 0:
        return a
    if sigma < EXACT_BELOW:
        return ndi.gaussian_filter1d(a, sigma, axis=axis, mode="nearest")
    w = int(round(np.sqrt(4.0 * sigma * sigma + 1.0)))
    if w % 2 == 0:
        w += 1
    out = a
    for _ in range(3):
        out = ndi.uniform_filter1d(out, w, axis=axis, mode="nearest")
    return out


def blur(a: np.ndarray, sigma) -> np.ndarray:
    """Separable version of `blur1d`; `sigma` is a scalar or one per axis."""
    sig = (sigma,) * a.ndim if np.isscalar(sigma) else tuple(sigma)
    out = a
    for ax, sg in enumerate(sig):
        out = blur1d(out, float(sg), ax)
    return out
