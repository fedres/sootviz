"""Finding the two families of triple-point trajectories.

A soot foil carries two families of trajectories, one running up-channel and
one down, and the lattice they form need not be symmetric: a sheared pattern
has the two families at different angles. The cell width is the TRANSVERSE
distance between adjacent trajectories of ONE family, so the families must be
separated before anything is measured.

They are separated with a bank of oriented line filters: each filter smooths
along its own direction and not across it, so a trajectory responds to the
filters aligned with it and to no others. The two modal angles of the bank are
the two families.

On a large record the pattern turns along the plate, so the angles are fitted
per block of a grid and the resulting fields are blended with a partition of
unity, which keeps the contours continuous across block boundaries.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy import ndimage as ndi

from .filters import blur, blur1d


def ridge_map(clean_img: np.ndarray, s: float, cell: float,
              polarity: str = "dark", fine_frac: float = 0.03,
              coarse_frac: float = 0.30) -> np.ndarray:
    """Trajectories -> positive. `polarity` is how tracks appear in the image.

    Numerical foils usually plot the maximum pressure, so the tracks are dark.
    Photographed foils may be either: the soot is scoured along a track (bright)
    or deposited there (dark). Use "auto" to let the bank decide.
    """
    fine = ndi.gaussian_filter(clean_img, fine_frac * cell * s)
    coarse = blur(clean_img, coarse_frac * cell * s)
    r = coarse - fine if polarity == "dark" else fine - coarse
    r = np.clip(r, 0, None)
    m = r.max()
    # float32 halves the filter bank, which is the largest thing in the
    # pipeline; the extra precision buys nothing on image data
    return (r / m if m > 0 else r).astype(np.float32)


def oriented(R: np.ndarray, deg: float, s: float, cell: float,
             length_frac: float = 0.30, width_frac: float = 0.03) -> np.ndarray:
    """Response of a line matched filter at `deg` (x right, y down-image)."""
    ny, nx = R.shape
    rot = ndi.rotate(R, -deg, reshape=True, order=1, mode="nearest")
    # across the line: exact, the shape is what selects a thin track.
    # along the line: wide, so the cheap approximation is used.
    sm = ndi.gaussian_filter1d(rot, width_frac * cell * s, axis=0,
                               mode="nearest")
    sm = blur1d(sm, length_frac * cell * s, axis=1)
    back = ndi.rotate(sm, deg, reshape=True, order=1, mode="nearest")
    cy, cx = back.shape[0] // 2, back.shape[1] // 2
    return back[cy - ny // 2:cy - ny // 2 + ny, cx - nx // 2:cx - nx // 2 + nx]


def bank(R: np.ndarray, s: float, cell: float, angles: np.ndarray,
         length_frac: float = 0.30, width_frac: float = 0.03) -> np.ndarray:
    return np.stack([oriented(R, a, s, cell, length_frac, width_frac)
                     for a in angles])


def family_angles(resp: np.ndarray, angles: np.ndarray,
                  region=None, min_abs: float = 10.0, max_abs: float = 65.0):
    """The two modal orientations, and a coherence score in [0, 1].

    The score is the fraction of strong track pixels that belong to the two
    modal directions (within one bank step). A pattern with two clean families
    scores high; noise scores low, and a score below about 0.25 means the
    two-family picture should not be trusted for that region.
    """
    best = angles[np.argmax(resp, axis=0)]
    strength = resp.max(axis=0)
    if region is not None:
        ys, ye, xs, xe = region
        best, strength = best[ys:ye, xs:xe], strength[ys:ye, xs:xe]
    sel = strength > np.percentile(strength, 70)
    if sel.sum() < 32:
        return None, None, 0.0
    counts = np.array([np.count_nonzero(best[sel] == a) for a in angles])
    pos = [i for i, a in enumerate(angles) if min_abs <= a <= max_abs]
    neg = [i for i, a in enumerate(angles) if -max_abs <= a <= -min_abs]
    if not pos or not neg or counts[pos].max() == 0 or counts[neg].max() == 0:
        return None, None, 0.0
    i1, i2 = pos[int(np.argmax(counts[pos]))], neg[int(np.argmax(counts[neg]))]
    step = 1
    near = lambda i: counts[max(0, i - step):i + step + 1].sum()
    score = (near(i1) + near(i2)) / max(1, counts.sum())
    return float(angles[i1]), float(angles[i2]), float(score)


@dataclass
class Block:
    xc: float
    yc: float
    a1: float
    a2: float
    score: float
    lam: float = float("nan")
    n: int = 0


def grid_blocks(shape, s, cell, block_frac=2.5, overlap=0.5, exclude=(),
                min_cells=4.0, min_height_frac=0.5, max_blocks=60):
    """Block centres covering the record, in physical units.

    A block supplies one pair of family angles, so it has to contain enough
    cells for the two modal directions to be separable: a block holding one or
    two cells returns the orientation of whichever few tracks it happens to
    catch. The requested size is therefore raised to hold at least `min_cells`
    cells and to span at least `min_height_frac` of the record height, and the
    count is capped so that a fine trial scale cannot produce thousands of
    blocks.
    """
    ny, nx = shape
    H_rec = ny / s
    W = max(block_frac * cell, min_cells * cell, min_height_frac * H_rec)
    W = min(W, nx / s)
    for _ in range(12):
        H = min(W, H_rec)
        step = max(W * (1 - overlap), 1e-12)
        xs = np.arange(W / 2, nx / s - W / 2 + 1e-9, step)
        if len(xs) == 0:
            xs = np.array([nx / s / 2])
        ys = np.arange(H / 2, H_rec - H / 2 + 1e-9, H) if H_rec > 1.5 * H \
            else np.array([H_rec / 2])
        if len(xs) * len(ys) <= max_blocks:
            break
        W *= 1.4
        W = min(W, nx / s)
    out = []
    for y in ys:
        for x in xs:
            if any(abs(x - e) < W / 2 for e in exclude):
                continue
            out.append((float(x), float(y), float(W), float(H)))
    return out


def blend(resp, angles, blocks, shape, s, band=10.0,
          hx_frac=0.7, hy_frac=1.2):
    """Blend the per-block family responses into two continuous fields."""
    ny, nx = shape
    X, Y = np.meshgrid(np.arange(nx) / s, np.arange(ny) / s)
    F1 = np.zeros(shape, np.float32); F2 = np.zeros(shape, np.float32)
    W = np.zeros(shape, np.float32)
    for b, (xc, yc, bw, bh) in blocks:
        if b is None:
            continue
        w = (np.clip(1 - np.abs(X - xc) / (hx_frac * bw), 0, None) *
             np.clip(1 - np.abs(Y - yc) / (hy_frac * bh), 0, None))
        if w.max() <= 0:
            continue
        for a0, F in ((b.a1, F1), (b.a2, F2)):
            sub = resp[(angles >= a0 - band) & (angles <= a0 + band)].max(axis=0)
            m = sub.max()
            F += w * (sub / m if m > 0 else sub)
        W += w
    good = W > 1e-6
    out = []
    for F in (F1, F2):
        G = np.where(good, F / np.where(good, W, 1), np.nan)
        out.append(G)
    return out[0], out[1], good


def local_normalise(F, s, cell, window_frac=2.0):
    """Divide by a running mean so one contour level suits the whole record."""
    G = np.nan_to_num(F)
    m = blur(G, window_frac * cell * s) + 1e-9
    return np.where(np.isnan(F), np.nan, G / m)
