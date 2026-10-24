"""Figures and data files, including the overlay the user must approve."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from .tracks import local_normalise

INK, GRID, SIG = "#222222", "#d9d9d9", "#b00020"
RED, BLUE = "#d81b3c", "#1f6feb"


def _mpl():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 8, "axes.linewidth": 0.5,
                         "figure.dpi": 150, "savefig.bbox": "tight",
                         "savefig.pad_inches": 0.03, "pdf.fonttype": 42})
    return plt


def fig_families(res, out, strips=None, level=1.25, dpi=170):
    """THE VERIFICATION FIGURE. Every contour should wrap exactly one
    trajectory; where contours merge lanes or wander off the structure, the
    measurement in that region is not to be trusted."""
    plt = _mpl()
    f = res.foil
    sa = res.s or f.s
    N1 = local_normalise(res.F1, sa, res.cell)
    N2 = local_normalise(res.F2, sa, res.cell)
    # the fields live on the analysis raster, which need not be the input
    # raster, so the axes come from the field itself
    X = np.linspace(0, f.width, N1.shape[1])
    Y = np.linspace(0, f.height, N1.shape[0])
    if strips is None:
        n = max(1, int(round((f.width / f.height) / 2.2)))
        edges = np.linspace(0, f.width, n + 1)
        strips = [(edges[i], edges[i + 1]) for i in range(n)]
    fig, axes = plt.subplots(len(strips), 1,
                             figsize=(13, 1.0 + 2.6 * len(strips)),
                             squeeze=False)
    for ax, (xa, xb) in zip(axes[:, 0], strips):
        ax.imshow(res.clean, cmap="gray", extent=f.extent, aspect="equal",
                  interpolation="lanczos", vmin=-2, vmax=2)
        ax.contour(X, Y, N1, levels=[level], colors=RED, linewidths=0.8)
        ax.contour(X, Y, N2, levels=[level], colors=BLUE, linewidths=0.8)
        ax.set_xlim(xa, xb); ax.set_ylim(f.height, 0)
        ax.set_ylabel(f"y [{f.unit}]"); ax.tick_params(labelsize=7)
    axes[-1, 0].set_xlabel(f"x [{f.unit}]")
    axes[0, 0].set_title(
        f"{f.name}: track families (red/blue) over the cleaned record — "
        f"each contour should wrap one trajectory", fontsize=9, loc="left")
    fig.savefig(Path(out).with_suffix(".png"), dpi=dpi)
    plt.close(fig)
    return Path(out).with_suffix(".png")


def fig_overview(res, out, level=1.25, dpi=200):
    """The whole record in one panel, contours blended across all blocks."""
    plt = _mpl()
    f = res.foil
    sa = res.s or f.s
    N1 = local_normalise(res.F1, sa, res.cell)
    N2 = local_normalise(res.F2, sa, res.cell)
    # the fields live on the analysis raster, which need not be the input
    # raster, so the axes come from the field itself
    X = np.linspace(0, f.width, N1.shape[1])
    Y = np.linspace(0, f.height, N1.shape[0])
    h = max(2.2, 13 * f.height / f.width + 0.9)
    fig, ax = plt.subplots(figsize=(13, h))
    ax.imshow(res.clean, cmap="gray", extent=f.extent, aspect="equal",
              interpolation="lanczos", vmin=-2, vmax=2)
    ax.contour(X, Y, N1, levels=[level], colors=RED, linewidths=0.6)
    ax.contour(X, Y, N2, levels=[level], colors=BLUE, linewidths=0.6)
    ax.set_xlabel(f"x [{f.unit}]"); ax.set_ylabel(f"y [{f.unit}]")
    ax.set_title(f"{f.name}: {f.width:g} x {f.height:g} {f.unit}, "
                 f"lambda = {res.median:.3g} {f.unit}", fontsize=9, loc="left")
    fig.savefig(Path(out).with_suffix(".png"), dpi=dpi)
    plt.close(fig)
    return Path(out).with_suffix(".png")


def _bands(reference):
    """Normalise whatever a caller passed as a reference into drawable bands."""
    if reference is None:
        return []
    if isinstance(reference, (int, float)):
        return [(float(reference), float(reference), "reference", BLUE)]
    if (isinstance(reference, (tuple, list)) and len(reference) == 2
            and all(isinstance(v, (int, float)) for v in reference)):
        return [(float(reference[0]), float(reference[1]), "reference", BLUE)]
    out = []
    for item in reference:
        (lo, hi), lab, col = item
        out.append((float(lo), float(hi), lab, col))
    return out


def fig_distribution(res, out, reference=None, dpi=170):
    """Histogram, and the block medians along the record."""
    plt = _mpl()
    f, st = res.foil, res.stats()
    lam = res.lam
    fig, axes = plt.subplots(1, 2, figsize=(8.4, 2.9))
    nb = max(12, min(40, int(np.sqrt(lam.size))))
    axes[0].hist(lam, bins=nb, color="#c9ccd1", edgecolor=INK, lw=0.4)
    axes[0].axvspan(st["q1"], st["q3"], color=SIG, alpha=0.10, lw=0)
    axes[0].axvline(st["median"], color=SIG, lw=1.3)
    axes[0].text(st["median"], axes[0].get_ylim()[1] * 0.94,
                 f"  median {st['median']:.3g}", color=SIG, fontsize=7.5)
    # `reference` is a single value, a (lo, hi) range, or a list of
    # ((lo, hi), label, colour) bands.
    for lo, hi, lab, col in _bands(reference):
        if hi > lo:
            axes[0].axvspan(lo, hi, color=col, alpha=0.12, lw=0)
        else:
            axes[0].axvline(lo, color=col, lw=1.1, ls=(0, (4, 2)))
        axes[0].text((lo + hi) / 2, axes[0].get_ylim()[1] * 0.99, lab,
                     color=col, fontsize=6.5, ha="center", va="top")
    axes[0].set_xlabel(f"cell width $\\lambda$ [{f.unit}]")
    axes[0].set_ylabel("spacings")
    axes[0].grid(True, axis="y", color=GRID, lw=0.4); axes[0].set_axisbelow(True)
    bl = [b for b, _ in res.blocks if b is not None and b.n > 2]
    if bl:
        axes[1].plot([b.xc for b in bl], [b.lam for b in bl], "o", ms=3.5,
                     color=SIG, mec="none", alpha=0.85)
        axes[1].axhline(st["median"], color=INK, lw=0.8, ls="--")
        axes[1].set_ylim(0, max(b.lam for b in bl) * 1.3)
    axes[1].set_xlabel(f"x along the record [{f.unit}]")
    axes[1].set_ylabel(f"block median $\\lambda$ [{f.unit}]")
    axes[1].grid(True, color=GRID, lw=0.4); axes[1].set_axisbelow(True)
    fig.suptitle(f"{f.name}", fontsize=8.5, y=1.02)
    fig.savefig(Path(out).with_suffix(".png"), dpi=dpi)
    plt.close(fig)
    return Path(out).with_suffix(".png")


def fig_scale(res, out, dpi=170):
    """The scale search: measured spacing against the scale assumed to measure
    it. A correct measurement is FLAT here -- the answer stops depending on the
    assumption. The chosen point is the knee, the coarsest scale that still
    resolves what it reports; to the left of it the filters are finer than the
    pattern and subdivide real spacings."""
    plt = _mpl()
    if not res.ladder:
        return None
    c = np.array([a for a, _ in res.ladder])
    l = np.array([b for _, b in res.ladder])
    u = res.foil.unit
    fig, ax = plt.subplots(figsize=(4.4, 3.0))
    lim = [min(c.min(), np.nanmin(l)) * 0.85, max(c.max(), np.nanmax(l)) * 1.1]
    ax.plot(lim, lim, color=GRID, lw=0.8, zorder=0)
    ax.plot(lim, [res.cfg.consistent * v for v in lim], color=GRID, lw=0.6,
            ls=":", zorder=0)
    ax.plot(c, l, "-o", color=INK, ms=3.4, lw=0.9, mfc="white", mew=0.9)
    ax.plot([res.cell], [res.median], "o", color=SIG, ms=7, mfc="none", mew=1.6,
            label=f"chosen: {res.median:.3g} {u}")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel(f"assumed scale c [{u}]")
    ax.set_ylabel(f"measured spacing [{u}]")
    ax.set_title("scale search: flat means the answer does not\n"
                 "depend on the assumption", fontsize=8, loc="left")
    ax.legend(frameon=False, fontsize=7, loc="lower right")
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    fig.savefig(Path(out).with_suffix(".png"), dpi=dpi)
    plt.close(fig)
    return Path(out).with_suffix(".png")


def write_report(res, outdir, reference=None, prefix=None):
    """All figures and data for one record. Returns the paths written."""
    outdir = Path(outdir); outdir.mkdir(parents=True, exist_ok=True)
    p = prefix or res.foil.name.replace(" ", "_")
    written = {
        "verification": fig_families(res, outdir / f"{p}_families"),
        "overview": fig_overview(res, outdir / f"{p}_overview"),
        "scale": fig_scale(res, outdir / f"{p}_scale"),
        "distribution": fig_distribution(res, outdir / f"{p}_distribution",
                                         reference=reference),
    }
    blocks = [dict(x=b.xc, y=b.yc, a1=b.a1, a2=b.a2, score=b.score,
                   lam=b.lam, n=b.n) for b, _ in res.blocks if b is not None]
    doc = dict(tool="SootViz-compute", name=res.foil.name,
               unit=res.foil.unit, width=res.foil.width,
               height=res.foil.height, px_per_unit=res.foil.s,
               px_per_unit_analysis=res.s,
               polarity=res.polarity, coherence=res.coherence,
               cell_scale_used=res.cell, stats=res.stats(),
               scale_search=[{"assumed": a, "measured": b}
                             for a, b in res.ladder],
               family_agreement_percent=res.family_agreement(),
               cells_across_height=res.cells_across(),
               scale_sensitivity_percent=res.scale_sensitivity(),
               family_medians=[float(np.median(res.lam1)) if len(res.lam1) else None,
                               float(np.median(res.lam2)) if len(res.lam2) else None],
               blocks=blocks)
    (outdir / f"{p}_result.json").write_text(json.dumps(doc, indent=1))
    np.savetxt(outdir / f"{p}_spacings.csv", res.lam, delimiter=",",
               header=f"cell width [{res.foil.unit}]", comments="")
    written = {k: v for k, v in written.items() if v is not None}
    written["json"] = outdir / f"{p}_result.json"
    written["csv"] = outdir / f"{p}_spacings.csv"
    return written
