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



def _fields(res):
    """Normalised family fields and their physical axes."""
    sa = res.s or res.foil.s
    N1 = local_normalise(res.F1, sa, res.cell)
    N2 = local_normalise(res.F2, sa, res.cell)
    X = np.linspace(0, res.foil.width, N1.shape[1])
    Y = np.linspace(0, res.foil.height, N1.shape[0])
    return N1, N2, X, Y


def fig_clean(res, out, dpi=170):
    """What the auto-contrast step did.

    The record as loaded, and after despeckling and local contrast
    equalisation over a window several cells wide. The window matters: equalise
    at the scale of a cell and the pattern being measured is removed with the
    illumination.
    """
    plt = _mpl()
    f = res.foil
    fig, ax = plt.subplots(2, 1, figsize=(12, 2.0 + 9.0 * f.height / f.width),
                           squeeze=False)
    ax[0, 0].imshow(f.image, cmap="gray", extent=f.extent, aspect="equal",
                    interpolation="lanczos")
    ax[0, 0].set_title(f"{f.name}: as loaded", fontsize=9, loc="left")
    ax[1, 0].imshow(res.clean, cmap="gray", extent=f.extent, aspect="equal",
                    interpolation="lanczos", vmin=-2, vmax=2)
    ax[1, 0].set_title(f"after despeckle and local contrast equalisation over "
                       f"{res.cfg.equalise_frac:g} cell widths "
                       f"({res.cfg.equalise_frac * res.cell:.3g} {f.unit})",
                       fontsize=9, loc="left")
    for a in ax[:, 0]:
        a.set_ylabel(f"y [{f.unit}]"); a.tick_params(labelsize=7)
    ax[1, 0].set_xlabel(f"x [{f.unit}]")
    fig.savefig(Path(out).with_suffix(".png"), dpi=dpi)
    plt.close(fig)
    return Path(out).with_suffix(".png")


def fig_blocks(res, out, level=1.25, max_panels=16, dpi=160):
    """The record as it was actually analysed: block by block.

    Each block supplies its OWN pair of family angles, because the pattern
    turns along a record. The panels show what each block saw and what it
    measured; the contours in the overview are these per-block fields blended
    with a partition of unity, which is why they stay continuous across the
    joins.
    """
    plt = _mpl()
    f, sa = res.foil, (res.s or res.foil.s)
    N1, N2, X, Y = _fields(res)
    rows = [(b, box) for b, box in res.blocks if b is not None]
    if not rows:
        return None
    note = ""
    if len(rows) > max_panels:                 # sample evenly, say so
        idx = np.linspace(0, len(rows) - 1, max_panels).round().astype(int)
        rows = [rows[i] for i in sorted(set(idx))]
        note = f" ({len(rows)} of {len(res.per_block())} blocks shown)"
    nc = min(4, len(rows))
    nr = int(np.ceil(len(rows) / nc))
    fig, axes = plt.subplots(nr, nc, figsize=(3.3 * nc, 3.0 * nr),
                             squeeze=False)
    for ax in axes.ravel():
        ax.axis("off")
    for ax, (b, (xc, yc, bw, bh)) in zip(axes.ravel(), rows):
        x0, x1 = max(0, xc - bw / 2), min(f.width, xc + bw / 2)
        y0, y1 = max(0, yc - bh / 2), min(f.height, yc + bh / 2)
        i0, i1 = int(y0 * sa), int(y1 * sa)
        j0, j1 = int(x0 * sa), int(x1 * sa)
        ax.axis("on")
        ax.imshow(res.clean[i0:i1, j0:j1], cmap="gray",
                  extent=(x0, x1, y1, y0), aspect="equal",
                  interpolation="lanczos", vmin=-2, vmax=2)
        ax.contour(X[j0:j1], Y[i0:i1], N1[i0:i1, j0:j1], levels=[level],
                   colors=RED, linewidths=0.7)
        ax.contour(X[j0:j1], Y[i0:i1], N2[i0:i1, j0:j1], levels=[level],
                   colors=BLUE, linewidths=0.7)
        lam = f"{b.lam:.3g}" if np.isfinite(b.lam) else "--"
        ax.set_title(f"x = {xc:.3g} {f.unit}\n"
                     f"$\\lambda$ = {lam}, n = {b.n}\n"
                     f"{b.a1:+.0f}$^\\circ$ / {b.a2:+.0f}$^\\circ$, "
                     f"coh {b.score:.2f}", fontsize=7.5, loc="left")
        ax.tick_params(labelsize=6)
    fig.suptitle(f"{f.name}: the blocks the angles were fitted in{note}",
                 fontsize=9.5, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.985))
    fig.savefig(Path(out).with_suffix(".png"), dpi=dpi)
    plt.close(fig)
    return Path(out).with_suffix(".png")


def fig_counting(res, out, n_cuts=9, dpi=180):
    """Where every number in the histogram came from.

    Spacings are read along vertical cuts through each family's field: a peak
    is one trajectory crossing the cut, and the gap between consecutive peaks
    is one cell width. Nothing else enters the statistics, so this figure can
    be checked by eye -- count the marks between two labels and the arithmetic
    is the whole method.
    """
    plt = _mpl()
    f, sa = res.foil, (res.s or res.foil.s)
    if not res.cuts1:
        return None
    # a window a few cells wide, in the middle of the record
    w = min(f.width, 6.0 * res.cell)
    xc = f.width / 2
    x0, x1 = max(0, xc - w / 2), min(f.width, xc + w / 2)
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 1.2 + 5.0 * f.height / w),
                             squeeze=False, sharey=True)
    for ax, cuts, col, lab in ((axes[0, 0], res.cuts1, RED, "family 1"),
                               (axes[0, 1], res.cuts2, BLUE, "family 2")):
        ax.imshow(res.clean, cmap="gray", extent=f.extent, aspect="auto",
                  interpolation="lanczos", vmin=-2, vmax=2)
        inside = [c for c in cuts if x0 <= c[0] <= x1 and len(c[1]) > 1]
        if inside:
            pick = [inside[i] for i in
                    np.linspace(0, len(inside) - 1, min(n_cuts, len(inside)))
                    .round().astype(int)]
            for k, (x, pks) in enumerate(dict.fromkeys(
                    (c[0], tuple(c[1])) for c in pick).keys()):
                pks = np.asarray(pks)
                ax.plot([x, x], [0, f.height], color=col, lw=0.4, alpha=0.35)
                ax.plot([x] * len(pks), pks, "o", color=col, ms=3.2,
                        mec="white", mew=0.5)
                if k == len(pick) // 2:          # label one cut's spacings
                    for a, bb in zip(pks[:-1], pks[1:]):
                        ax.annotate(f"{bb - a:.3g}", (x, (a + bb) / 2),
                                    xytext=(5, 0), textcoords="offset points",
                                    fontsize=6.5, color=col, va="center")
                        ax.annotate("", (x, bb), (x, a), arrowprops=dict(
                            arrowstyle="<->", color=col, lw=0.7))
        ax.set_xlim(x0, x1); ax.set_ylim(f.height, 0)
        ax.set_xlabel(f"x [{f.unit}]")
        ax.set_title(f"{lab}: {len(cuts)} cuts over the record, "
                     f"{len(res.lam1 if col == RED else res.lam2)} spacings",
                     fontsize=8.5, loc="left")
    axes[0, 0].set_ylabel(f"y [{f.unit}]")
    fig.suptitle(f"{f.name}: how the spacings were counted "
                 f"(every gap between consecutive marks is one measurement)",
                 fontsize=9.5, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(Path(out).with_suffix(".png"), dpi=dpi)
    plt.close(fig)
    return Path(out).with_suffix(".png")


def fig_distribution(res, out, reference=None, dpi=170):
    """The result is this distribution, not the single number taken from it.

    Three panels: every spacing measured; the same split into the two
    families, which are independent measurements of the same quantity, so
    their overlap is the internal check; and the block medians along the
    record, which show whether the cell size drifts downstream.
    """
    plt = _mpl()
    f, st = res.foil, res.stats()
    lam, u = res.lam, res.foil.unit
    fig, axes = plt.subplots(1, 3, figsize=(12.6, 3.1))

    nb = max(14, min(48, int(np.sqrt(lam.size))))
    edges = np.histogram_bin_edges(lam, bins=nb)
    axes[0].hist(lam, bins=edges, color="#c9ccd1", edgecolor=INK, lw=0.4)
    axes[0].axvspan(st["q1"], st["q3"], color=SIG, alpha=0.10, lw=0)
    axes[0].axvline(st["median"], color=SIG, lw=1.4)
    axes[0].text(st["median"], axes[0].get_ylim()[1] * 0.95,
                 f"  median {st['median']:.3g}", color=SIG, fontsize=7.5)
    for lo, hi, lab, col in _bands(reference):
        if hi > lo:
            axes[0].axvspan(lo, hi, color=col, alpha=0.12, lw=0)
        else:
            axes[0].axvline(lo, color=col, lw=1.1, ls=(0, (4, 2)))
        axes[0].text((lo + hi) / 2, axes[0].get_ylim()[1] * 0.99, lab,
                     color=col, fontsize=6.5, ha="center", va="top")
    axes[0].set_title(f"all {st['n']} spacings   IQR "
                      f"{st['q1']:.3g}-{st['q3']:.3g} {u}",
                      fontsize=8.5, loc="left")
    axes[0].set_xlabel(f"cell width $\\lambda$ [{u}]")
    axes[0].set_ylabel("spacings")

    for a, col, lab in ((res.lam1, RED, "family 1"), (res.lam2, BLUE, "family 2")):
        if a is not None and len(a):
            axes[1].hist(a, bins=edges, histtype="step", lw=1.2, color=col,
                         label=f"{lab}: {np.median(a):.3g} {u} (n={len(a)})")
    axes[1].legend(frameon=False, fontsize=7)
    axes[1].set_title(f"the two families, measured separately "
                      f"— they differ by {res.family_agreement():.0f} %",
                      fontsize=8.5, loc="left")
    axes[1].set_xlabel(f"cell width $\\lambda$ [{u}]")

    bl = [b for b, _ in res.blocks if b is not None and b.n > 2]
    if bl:
        n = np.array([b.n for b in bl], float)
        axes[2].scatter([b.xc for b in bl], [b.lam for b in bl],
                        s=8 + 40 * n / n.max(), color=SIG, alpha=0.8,
                        edgecolors="none")
        axes[2].axhline(st["median"], color=INK, lw=0.8, ls="--")
        axes[2].set_ylim(0, max(b.lam for b in bl) * 1.35)
        axes[2].set_title(f"{len(bl)} blocks along the record "
                          f"(marker area = spacings counted)",
                          fontsize=8.5, loc="left")
    axes[2].set_xlabel(f"x along the record [{u}]")
    axes[2].set_ylabel(f"block median $\\lambda$ [{u}]")

    for a in axes:
        a.grid(True, color=GRID, lw=0.4); a.set_axisbelow(True)
        for sp in ("top", "right"):
            a.spines[sp].set_visible(False)
    fig.suptitle(f"{f.name}", fontsize=8.5, y=1.03, x=0.005, ha="left")
    fig.tight_layout()
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


def _json_safe(o):
    """NaN and infinity are not JSON.

    Python writes them as bare NaN/Infinity, which every strict parser
    rejects -- JavaScript's JSON.parse among them, so a single missing
    statistic would discard the whole result. They become null instead, which
    is what "not measured" means here: a record with one empty family has no
    family agreement, and a run given its scale has no scale search to be
    sensitive to.
    """
    if isinstance(o, dict):
        return {k: _json_safe(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_json_safe(v) for v in o]
    if isinstance(o, (float, np.floating)):
        return float(o) if np.isfinite(o) else None
    if isinstance(o, (int, np.integer)):
        return int(o)
    return o


def write_report(res, outdir, reference=None, prefix=None):
    """All figures and data for one record. Returns the paths written."""
    outdir = Path(outdir); outdir.mkdir(parents=True, exist_ok=True)
    p = prefix or res.foil.name.replace(" ", "_")
    written = {
        "verification": fig_families(res, outdir / f"{p}_families"),
        "overview": fig_overview(res, outdir / f"{p}_overview"),
        "clean": fig_clean(res, outdir / f"{p}_clean"),
        "blocks": fig_blocks(res, outdir / f"{p}_blocks"),
        "counting": fig_counting(res, outdir / f"{p}_counting"),
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
    (outdir / f"{p}_result.json").write_text(
        json.dumps(_json_safe(doc), indent=1, allow_nan=False))
    u = res.foil.unit
    fam = np.concatenate([np.ones(len(res.lam1)), 2 * np.ones(len(res.lam2))])
    tab = np.column_stack([res.lam, res.lam_x, fam]) \
        if len(res.lam_x) == len(res.lam) else res.lam[:, None]
    np.savetxt(outdir / f"{p}_spacings.csv", tab, delimiter=",", fmt="%.6g",
               header=f"cell_width[{u}],x[{u}],family", comments="")
    with open(outdir / f"{p}_blocks.csv", "w") as fh:
        fh.write(f"x[{u}],y[{u}],width[{u}],height[{u}],angle1[deg],"
                 f"angle2[deg],coherence,lambda[{u}],spacings\n")
        for b in res.per_block():
            fh.write(f"{b['x']:.6g},{b['y']:.6g},{b['width']:.6g},"
                     f"{b['height']:.6g},{b['a1']:.6g},{b['a2']:.6g},"
                     f"{b['coherence']:.4g},{b['lam']:.6g},{b['n']}\n")
    written = {k: v for k, v in written.items() if v is not None}
    written["json"] = outdir / f"{p}_result.json"
    written["csv"] = outdir / f"{p}_spacings.csv"
    written["blocks_csv"] = outdir / f"{p}_blocks.csv"
    return written
