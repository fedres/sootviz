"""Measuring the cell width from the separated track families."""
from __future__ import annotations

from dataclasses import dataclass, field, asdict

import numpy as np
from scipy import ndimage as ndi
from scipy.signal import find_peaks
from scipy.stats import skew

from .clean import clean as clean_image
from .io import Foil
from .tracks import (Block, bank, blend, family_angles, grid_blocks,
                     local_normalise, ridge_map)


@dataclass
class Config:
    """Every length is a fraction of the cell width, so the pipeline is
    scale-free: the same numbers work on a 30 mm numerical foil and a 5 m
    photographed plate."""
    cell: float | None = None        # expected cell width; None -> auto
    polarity: str = "auto"           # "dark", "bright" or "auto"
    angles: tuple = (-70, 70, 5)     # bank: min, max, step in degrees
    block_frac: float = 2.5          # grid block size, in cells
    grid: bool = True                # False -> one global pair of angles
    exclude_x: tuple = ()            # physical x to skip (seams, damage)
    # cleaning
    despeckle_frac: float = 0.03
    equalise_frac: float = 2.5
    denoise_frac: float = 0.05
    # ridge map and filter bank
    fine_frac: float = 0.024
    coarse_frac: float = 0.20
    length_frac: float = 0.24
    width_frac: float = 0.017
    # sampling along the cuts
    cut_frac: float = 0.10           # spacing between vertical cuts
    smooth_frac: float = 0.02        # smoothing along a cut
    min_sep_frac: float = 0.22       # minimum separation of two trajectories
    prominence: float = 0.12
    # scale search
    ladder: float = 0.8             # step of the descending scale ladder
    consistent: float = 0.85        # lambda_hat >= this x c counts as resolved
    # The analysis runs at a canonical resolution of this many pixels per cell.
    # Filter sizes are fractions of a cell, so without this the cost grows with
    # the raster resolution and with the trial scale -- a coarse trial scale on
    # a high-resolution scan needs Gaussian kernels hundreds of pixels wide for
    # no gain in information. Decimating first makes the cost, and the result,
    # independent of how finely the record happens to be sampled.
    #
    # The value is set from below by the thinnest filter: the cross-track width
    # is 0.020 of a cell, so 96 px per cell puts it at about 2 px. Coarser than
    # that and the cross-track filter goes sub-pixel, stops rejecting noise,
    # and the scale search runs away. (2 px across, 28 px along, is also the
    # filter geometry of the hand measurement this pipeline reproduces.)
    px_per_cell: float = 96.0


@dataclass
class Result:
    foil: Foil
    cfg: Config
    cell: float
    lam: np.ndarray
    blocks: list
    polarity: str
    coherence: float
    # pixels per unit of the ANALYSIS raster, which may be a decimated copy of
    # the input: the fields and the cleaned image are at this resolution.
    s: float = 0.0
    clean: np.ndarray = field(repr=False, default=None)
    F1: np.ndarray = field(repr=False, default=None)
    F2: np.ndarray = field(repr=False, default=None)
    lam1: np.ndarray = field(repr=False, default=None)
    lam2: np.ndarray = field(repr=False, default=None)
    # (assumed scale, measured spacing) over the search, for the scale panel
    ladder: list = field(repr=False, default_factory=list)

    # ---- statistics -------------------------------------------------------
    @property
    def median(self):
        return float(np.median(self.lam)) if len(self.lam) else float("nan")

    def stats(self):
        a = self.lam
        if not len(a):
            return {}
        q1, q3 = np.percentile(a, [25, 75])
        return dict(n=int(a.size), median=float(np.median(a)),
                    mean=float(a.mean()), sd=float(a.std(ddof=1)),
                    q1=float(q1), q3=float(q3),
                    p5=float(np.percentile(a, 5)),
                    p95=float(np.percentile(a, 95)), max=float(a.max()))

    def family_agreement(self):
        """Per-cent difference between the two families' medians. The two are
        independent measurements of the same quantity; disagreement above
        about 15 % means the geometry has not been read correctly."""
        if self.lam1 is None or self.lam2 is None or not len(self.lam1) \
                or not len(self.lam2):
            return float("nan")
        m1, m2 = np.median(self.lam1), np.median(self.lam2)
        return float(100 * abs(m1 - m2) / np.median(self.lam))

    def scale_sensitivity(self):
        """How much the answer moves with the assumed scale, in per cent.

        Taken over the chosen scale and the two coarser steps of the search.
        Near zero means the measurement sits on a plateau and the cell width is
        a sharp quantity for this record. A large value is not necessarily an
        error: an irregular front carries structure over a range of scales, and
        then the number quoted should be accompanied by the scale it was taken
        at."""
        if not self.ladder:
            return float("nan")
        cs = [c for c, _ in self.ladder]
        try:
            i = cs.index(self.cell)
        except ValueError:
            return float("nan")
        near = [l for _, l in self.ladder[max(0, i - 2):i + 1]
                if np.isfinite(l)]
        if len(near) < 2:
            return float("nan")
        return float(100 * (max(near) - min(near)) / np.median(near))

    def cells_across(self):
        return self.foil.height / self.median

    def summary(self):
        st = self.stats()
        u = self.foil.unit
        if not st:
            return "no spacings measured"
        return (f"{self.foil.name}: lambda = {st['median']:.3g} {u} "
                f"[IQR {st['q1']:.3g}-{st['q3']:.3g}, n={st['n']}]\n"
                f"  {self.cells_across():.1f} cells across the "
                f"{self.foil.height:g} {u} height\n"
                f"  families agree to {self.family_agreement():.0f} %, "
                f"coherence {self.coherence:.2f}, polarity '{self.polarity}', "
                f"{len(self.blocks)} block(s)"
                + ("" if not np.isfinite(self.scale_sensitivity()) else
                   f"\n  scale sensitivity {self.scale_sensitivity():.0f} % "
                   f"(near zero = the answer does not depend on the "
                   f"assumed scale)"))


def _spacings(F, foil: Foil, cfg: Config, cell: float, region=None,
              with_x: bool = False):
    s = foil.s
    ny, nx = F.shape
    ys, ye, xs, xe = region if region else (0, ny, 0, nx)
    step = max(1, foil.px(cfg.cut_frac * cell))
    sm = cfg.smooth_frac * cell * s
    dist = max(2, foil.px(cfg.min_sep_frac * cell))
    out, xo = [], []
    for xp in range(xs, xe, step):
        col = F[ys:ye, xp]
        if not np.isfinite(col).all() or col.max() <= 0:
            continue
        col = ndi.gaussian_filter1d(col, sm)
        pk, _ = find_peaks(col, prominence=cfg.prominence * col.max(),
                           distance=dist)
        if len(pk) > 1:
            d = list(np.diff(np.sort(pk)) / s)
            out += d
            xo += [xp / s] * len(d)
    if with_x:
        return np.asarray(out), np.asarray(xo)
    return np.asarray(out)


def _decimate(foil: Foil, cell: float, px_per_cell: float) -> Foil:
    """Resample so that one cell spans about `px_per_cell` pixels."""
    k = int(foil.px(cell) / px_per_cell)
    if k < 2:
        return foil
    img = ndi.uniform_filter(foil.image, k)[::k, ::k]   # mean, then take
    return Foil(img, foil.s / k, foil.width, foil.height, foil.unit, foil.name)


def _once(foil: Foil, cfg: Config, cell: float, polarity: str,
          quick: bool = False):
    foil = _decimate(foil, cell, cfg.px_per_cell)
    s = foil.s
    cl = clean_image(foil.image, s, cell, cfg.despeckle_frac,
                     cfg.equalise_frac, cfg.denoise_frac)
    R = ridge_map(cl, s, cell, polarity, cfg.fine_frac, cfg.coarse_frac)
    lo, hi, st = cfg.angles
    angles = np.arange(lo, hi + 1e-9, st)
    resp = bank(R, s, cell, angles, cfg.length_frac, cfg.width_frac)
    ny, nx = foil.image.shape

    # The grid fits the LOCAL family angles, because the pattern turns along
    # a long record. The spacings are then measured on the blended fields over
    # full-height cuts: measuring inside a block would bias the result low,
    # since peak prominence is judged against each cut's own maximum and a
    # short cut makes weak maxima look prominent.
    # The scale search only needs the median spacing, so it runs with one
    # global pair of angles: the grid and the per-block trend are for the
    # final measurement.
    blocks = []
    if cfg.grid and not quick:
        boxes = grid_blocks((ny, nx), s, cell, cfg.block_frac,
                            exclude=cfg.exclude_x)
    else:
        boxes = [(nx / s / 2, ny / s / 2, nx / s, ny / s)]
    for xc, yc, bw, bh in boxes:
        reg = (max(0, foil.px(yc - bh / 2)), min(ny, foil.px(yc + bh / 2)),
               max(0, foil.px(xc - bw / 2)), min(nx, foil.px(xc + bw / 2)))
        a1, a2, score = family_angles(resp, angles, region=reg)
        blocks.append((Block(xc, yc, a1, a2, score) if a1 is not None else None,
                       (xc, yc, bw, bh)))

    F1, F2, good = blend(resp, angles, blocks, (ny, nx), s, band=st)
    lam1, x1 = _spacings(np.nan_to_num(F1), foil, cfg, cell, with_x=True)
    lam2, x2 = _spacings(np.nan_to_num(F2), foil, cfg, cell, with_x=True)
    lam1, lam2 = list(lam1), list(lam2)

    if not quick:
        # the trend along the record, binned from the cuts already taken
        xall = np.concatenate([x1, x2]) if len(x1) or len(x2) else np.zeros(0)
        lall = np.asarray(lam1 + lam2)
        for b, (xc, yc, bw, bh) in blocks:
            if b is None:
                continue
            m = np.abs(xall - xc) <= bw / 2 if len(xall) else np.zeros(0, bool)
            b.lam = float(np.median(lall[m])) if m.any() else float("nan")
            b.n = int(m.sum())

    lam = np.asarray(lam1 + lam2)
    goodb = [b for b, _ in blocks if b is not None]
    coherence = float(np.median([b.score for b in goodb])) if goodb else 0.0
    return dict(clean=cl, s=s, lam=lam, lam1=np.asarray(lam1), lam2=np.asarray(lam2),
                blocks=blocks, coherence=coherence, F1=F1, F2=F2)


def _polarity(foil: Foil, cfg: Config, cell: float) -> tuple:
    """Decide whether the trajectories are dark or bright in this image.

    Trajectories are thin; what lies between them is the wide interior of a
    cell. The ridge map built with the correct sign is therefore concentrated
    on a small part of the area -- a positively skewed, heavy-tailed
    distribution of values -- while the wrong sign spreads its response over
    the interiors. The skewness of the ridge map separates the two cleanly
    (typically 2.2 against 0.8) and costs two filter passes, no bank.

    Numerical foils plotting maximum pressure come out "dark"; a photographed
    plate where the soot is scoured along a track comes out "bright".
    """
    d = _decimate(foil, cell, cfg.px_per_cell)
    cl = clean_image(d.image, d.s, cell, cfg.despeckle_frac,
                     cfg.equalise_frac, cfg.denoise_frac)
    sc = {}
    for pol in ("dark", "bright"):
        R = ridge_map(cl, d.s, cell, pol, cfg.fine_frac, cfg.coarse_frac)
        sc[pol] = float(skew(R.ravel()))
    return max(sc, key=sc.get), sc


def _solve(foil: Foil, cfg: Config, pol: str, verbose=False):
    """Find the coarsest scale at which the measurement resolves what it reports.

    Every length in the pipeline is a fraction of the cell width, which is the
    quantity being measured, so the scale has to be chosen. Measured against
    the assumed scale c, lambda_hat(c) behaves in two regimes: once the filters
    are fine enough to separate adjacent trajectories lambda_hat is flat -- the
    answer stops depending on the assumption, which is the signature of being
    right -- and below that a filter finer than the pattern starts subdividing
    real spacings, so lambda_hat follows c downwards. Coming down from coarse,
    the knee between the two is the first c whose measured spacing is
    consistent with the scale used to measure it.

    Two simpler rules were tried and are both wrong. Fixed-point iteration
    stops wherever its own tolerance happens to bite. Root-finding on
    lambda_hat(c) - c looks principled but on a regular lattice there is no
    root near the answer -- lambda_hat stays flat at 2.5 while c falls past it
    -- and the only root lies deep in the subdividing regime, which returned
    1.5 mm for a lattice built with 2.5 mm cells.
    """
    rows, best = [], None
    c, stop = foil.height / 3.0, foil.height / 20.0
    while c >= stop:
        out = _once(foil, cfg, c, pol, quick=True)
        lam = float(np.median(out["lam"])) if len(out["lam"]) else float("nan")
        rows.append((float(c), lam))
        if verbose:
            flag = "resolved" if np.isfinite(lam) and lam >= cfg.consistent * c \
                else "subdividing" if np.isfinite(lam) else "nothing found"
            print(f"    c = {c:9.4g} -> lambda = {lam:9.4g}  "
                  f"({lam / c:4.2f} of c, {flag})")
        if np.isfinite(lam) and lam >= cfg.consistent * c:
            best = (lam, float(c))
            break
        c *= cfg.ladder
    if best is None:
        fin = [(l / cc, l, cc) for cc, l in rows if np.isfinite(l)]
        if not fin:
            raise RuntimeError("no trajectories found at any scale: check the "
                               "physical extent, or set Config(cell=...)")
        _, l, cc = max(fin)
        best = (l, cc)
    return best[0], best[1], rows


def measure(foil: Foil, cfg: Config | None = None, verbose: bool = False) -> Result:
    """Measure the cell width, choosing the scale and the polarity if needed."""
    cfg = cfg or Config()

    pol = cfg.polarity
    if pol == "auto":
        pol, sc = _polarity(foil, cfg, cfg.cell or foil.height / 8.0)
        if verbose:
            print("  polarity: " + ", ".join(f"{k} skew {v:.2f}"
                                             for k, v in sc.items())
                  + f" -> '{pol}'")

    ladder = []
    if cfg.cell:
        out = _once(foil, cfg, cfg.cell, pol)
        cell = cfg.cell
    else:
        _, cell, ladder = _solve(foil, cfg, pol, verbose)
        out = _once(foil, cfg, cell, pol)      # final pass, with the grid
    if not len(out["lam"]):
        raise RuntimeError("no trajectories found: check the physical extent, "
                           "or set Config(cell=...) to the expected cell width")
    if verbose:
        print(f"    -> lambda = {np.median(out['lam']):.4g} {foil.unit}, "
              f"coherence {out['coherence']:.2f}")
    return Result(foil=foil, cfg=cfg, cell=cell, lam=out["lam"], s=out["s"],
                  blocks=out["blocks"], polarity=pol,
                  coherence=out["coherence"], clean=out["clean"],
                  F1=out["F1"], F2=out["F2"],
                  lam1=out["lam1"], lam2=out["lam2"], ladder=ladder)
