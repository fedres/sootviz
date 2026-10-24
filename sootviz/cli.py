"""Command line for SootViz-compute."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from . import __version__
from .io import list_pdf_rasters, load_foil, synth_foil
from .measure import Config, measure
from .report import write_report

BANNER = f"SootViz-compute {__version__} — Karthik Thyagarajan, Strawhat CFD"


def _add_common(p):
    p.add_argument("--width", type=float, required=True,
                   help="physical width of the image")
    p.add_argument("--height", type=float, required=True,
                   help="physical height of the image (the channel height)")
    p.add_argument("--unit", default="mm")
    p.add_argument("--cell", type=float, default=None,
                   help="expected cell width; omitted, it is found by iteration")
    p.add_argument("--polarity", choices=("auto", "dark", "bright"),
                   default="auto", help="how the trajectories appear")
    p.add_argument("--no-grid", action="store_true",
                   help="fit one pair of family angles for the whole record")
    p.add_argument("--block", type=float, default=2.5,
                   help="grid block size, in cell widths (default 2.5)")
    p.add_argument("--equalise", type=float, default=2.5,
                   help="contrast-equalisation window, in cell widths "
                        "(default 2.5; a window near one cell erases the "
                        "pattern, so raise it rather than lower it)")
    p.add_argument("--exclude-x", type=float, nargs="*", default=(),
                   help="x positions to skip (plate seams, damage)")
    p.add_argument("--pdf-page", type=int, default=None)
    p.add_argument("--pdf-index", type=int, default=0)
    p.add_argument("--invert", action="store_true")
    p.add_argument("-o", "--outdir", default="sootviz_out")


def cmd_measure(a):
    foil = load_foil(a.image, a.width, a.height, a.unit,
                     pdf_page=a.pdf_page, pdf_index=a.pdf_index,
                     invert=a.invert)
    cfg = Config(cell=a.cell, polarity=a.polarity, grid=not a.no_grid,
                 block_frac=a.block, exclude_x=tuple(a.exclude_x),
                 equalise_frac=a.equalise)
    print(BANNER); print(f"  {foil}")
    res = measure(foil, cfg, verbose=True)
    print(); print(res.summary())
    out = write_report(res, a.outdir)
    print("\nwritten:")
    for k, v in out.items():
        print(f"  {k:13s} {v}")
    print("\nCHECK THE VERIFICATION FIGURE BEFORE USING THE NUMBER.")
    print("  Every contour should wrap exactly one trajectory. Where contours")
    print("  merge neighbouring lanes or wander off the structure, the")
    print("  measurement there is not reliable: re-run with --cell set near")
    print("  the size you can see, or restrict the record with --exclude-x.")
    if res.family_agreement() > 15:
        print(f"\n  WARNING: the two families disagree by "
              f"{res.family_agreement():.0f} % (above 15 %), so the geometry "
              f"may have been read wrongly.")
    import math
    if (res.scale_sensitivity() or 0) > 25 and \
            math.isfinite(res.scale_sensitivity()):
        print(f"\n  NOTE: the answer moves {res.scale_sensitivity():.0f} % "
              f"across neighbouring scales, so this record has no single cell "
              f"width. Quote the distribution and the scale, not the median "
              f"alone -- see the scale panel.")
    if res.coherence < 0.25:
        print(f"\n  WARNING: orientation coherence {res.coherence:.2f} is low; "
              f"the record may not show two clean families.")
    return 0


def cmd_synth(a):
    foil = synth_foil(width=a.width, height=a.height, lam=a.lam, s=a.s,
                      skew=a.skew, jitter=a.jitter, noise=a.noise,
                      contrast=a.contrast, vignette=a.vignette, unit=a.unit)
    from PIL import Image
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    Image.fromarray(foil.image.astype("uint8")).save(a.out)
    print(f"{BANNER}\n  wrote {a.out}: {foil}, true lambda = {a.lam} {a.unit}")
    return 0


def cmd_benchmark(a):
    from .benchmark import run_published, run_synthetic
    print(BANNER)
    print("\nsynthetic foils (cell width known by construction):")
    ok, rows = run_synthetic(outdir=a.out)
    allrows = {"synthetic": rows}
    if a.pdf:
        print(f"\npublished foils from {a.pdf}:")
        ok2, rows2 = run_published(a.pdf, outdir=a.out)
        allrows["published"] = rows2
        ok = ok and ok2
    if a.out:
        import json
        p = Path(a.out); p.mkdir(parents=True, exist_ok=True)
        (p / "benchmarks.json").write_text(json.dumps(allrows, indent=1))
        print(f"\ngallery written to {p}")
    print("\nall benchmarks passed" if ok else "\nSOME BENCHMARKS FAILED")
    return 0 if ok else 1


def cmd_rasters(a):
    for i, size, aspect in list_pdf_rasters(a.pdf, a.page):
        print(f"  index {i}: {size[0]}x{size[1]} px, aspect {aspect:.2f}")
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(prog="sootviz", description=BANNER)
    ap.add_argument("--version", action="version", version=BANNER)
    sub = ap.add_subparsers(dest="cmd", required=True)

    m = sub.add_parser("measure", help="measure a soot foil")
    m.add_argument("image", help="image file, or a PDF page with --pdf-page")
    _add_common(m); m.set_defaults(func=cmd_measure)

    s = sub.add_parser("synth", help="write a synthetic foil of known cell width")
    s.add_argument("out"); s.add_argument("--lam", type=float, default=2.5)
    s.add_argument("--width", type=float, default=30.0)
    s.add_argument("--height", type=float, default=20.0)
    s.add_argument("--s", type=float, default=25.0)
    s.add_argument("--skew", type=float, default=0.0)
    s.add_argument("--jitter", type=float, default=0.12)
    s.add_argument("--noise", type=float, default=0.25)
    s.add_argument("--contrast", type=float, default=0.6)
    s.add_argument("--vignette", type=float, default=0.5)
    s.add_argument("--unit", default="mm"); s.set_defaults(func=cmd_synth)

    b = sub.add_parser("benchmark", help="run the benchmarks")
    b.add_argument("--pdf", default=None,
                   help="a published paper to measure (your own copy)")
    b.add_argument("--out", default=None,
                   help="write the figures and a summary JSON here, for the "
                        "benchmark gallery on the web page")
    b.set_defaults(func=cmd_benchmark)

    r = sub.add_parser("rasters", help="list the rasters on a PDF page")
    r.add_argument("pdf"); r.add_argument("--page", type=int, default=None)
    r.set_defaults(func=cmd_rasters)

    a = ap.parse_args(argv)
    return a.func(a)


if __name__ == "__main__":
    sys.exit(main())
