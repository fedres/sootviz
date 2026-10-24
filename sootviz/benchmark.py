"""Benchmarks: does the pipeline return the right number?

Three levels of evidence, in increasing order of realism:

  synthetic   a lattice whose cell width is known exactly by construction.
              The only test where truth is known; run it first after any
              change to the pipeline.
  published   soot foils from the literature, measured from the published
              figure. The expected values are what those authors report, so
              agreement shows the method is calibrated against the community's
              own numbers.
  own         records whose measurement has been published with this tool.
"""
from __future__ import annotations

import numpy as np

from .io import load_foil, synth_foil
from .measure import Config, measure

# Each case carries its own tolerance, because the cases are not equally
# measurable and a single pass mark would hide that. The scale search is tuned
# for real records, where the risk is drifting from the cell onto the GROUPING
# of cells, so it prefers the finer end of the plateau; on a heavily degraded
# image the plateau is short and that preference costs about 15 %. Reported,
# not hidden: a real record in that state also reports a high scale
# sensitivity, which is the signal to quote a distribution and not a number.
SYNTHETIC = [
    dict(name="regular", lam=2.5, aspect=1.8, skew=0.0, jitter=0.05,
         noise=0.15, tol=10.0),
    dict(name="sheared", lam=2.5, aspect=1.8, skew=12.0, jitter=0.08,
         noise=0.20, tol=10.0),
    dict(name="irregular", lam=2.0, aspect=1.6, skew=6.0, jitter=0.20,
         noise=0.35, tol=10.0),
    dict(name="noisy, low contrast", lam=3.0, aspect=1.8, skew=0.0,
         jitter=0.12, noise=0.60, contrast=0.35, vignette=0.7, tol=20.0),
]

# Published records. Paths are not bundled: the PDFs belong to their
# publishers. Point --pdf at your own copy.
#
# What the authors actually report matters here, so it is quoted rather than
# paraphrased. Taileb, Melguizo-Gavilanes & Chinnayya, Combust. Flame 218
# (2020), p. 17:
#
#   "The cell size was numerically estimated to lie in the range of
#    0.9 mm <= lam_num <= 1.9 mm whereas the cell size reported
#    experimentally is in the range of 1.4 mm <= lam_exp <= 2.1 mm"
#
# and, on the same page, "The smaller cell size predicted by detailed
# chemistry", which fixes the ordering of the three panels. Single-step was
# tuned to the experiment (p. 10: "by varying As and Ea/Ru to match the
# numerically determined cell size and steady detonation velocities with those
# reported experimentally"), so it sits at the top of the range.
#
# The authors give a range over the three mechanisms, not a value per panel.
# The honest check is therefore: every panel inside the reported range, and the
# ordering they state reproduced. Per-panel numbers are reported, not asserted.
PUBLISHED_RANGE = dict(
    source="Taileb, Melguizo-Gavilanes & Chinnayya, Combust. Flame 218 (2020) 247",
    numerical=(0.9, 1.9), experimental=(1.4, 2.1), unit="mm", page=17,
)
PUBLISHED = [
    dict(key="taileb2020-single-step", page=19, index=0, width=30.0,
         height=20.0, unit="mm", rank=0,
         cite="Fig. 5 top; single-step chemistry. Tuned by the authors to the "
              "experimental cell size, so expected at the top of the range"),
    dict(key="taileb2020-three-step", page=19, index=1, width=30.0,
         height=20.0, unit="mm", rank=1,
         cite="Fig. 5 centre; three-step chain-branching chemistry"),
    dict(key="taileb2020-detailed", page=19, index=2, width=30.0,
         height=20.0, unit="mm", rank=2,
         cite="Fig. 5 bottom; detailed chemistry. The authors note this gives "
              "the smaller cell size, so expected at the bottom of the range"),
]


def _write(res, outdir, prefix, reference=None):
    from .report import write_report
    return write_report(res, outdir, reference=reference, prefix=prefix)


def run_synthetic(verbose=True, seed=0, outdir=None):
    """Measure synthetic foils of known cell width. Returns (passed, rows)."""
    rows = []
    for case in SYNTHETIC:
        kw = dict(case); name = kw.pop("name"); tol = kw.pop("tol", 10.0)
        truth = kw["lam"]
        foil = synth_foil(seed=seed, **kw)
        foil.name = f"synthetic/{name}"
        res = measure(foil)
        if outdir:
            _write(res, outdir, f"synthetic-{name.replace(', ', '-').replace(' ', '-')}", truth)
        err = 100 * (res.median - truth) / truth
        ok = abs(err) <= tol
        rows.append(dict(name=name, truth=truth, measured=res.median,
                         error_percent=err, agreement=res.family_agreement(),
                         coherence=res.coherence, passed=ok))
        if verbose:
            print(f"  {name:22s} true {truth:5.2f}  measured "
                  f"{res.median:5.2f}  ({err:+5.1f} %)  "
                  f"families {res.family_agreement():4.0f} %  "
                  f"coherence {res.coherence:.2f}  "
                  f"{'ok' if ok else 'FAIL'} (tol {tol:.0f} %)")
    return all(r["passed"] for r in rows), rows


def run_published(pdf, keys=None, verbose=True, tol=0.25, outdir=None):
    """Measure published foils from a PDF the user supplies.

    The authors report a range over mechanisms rather than a value per panel,
    so two things are checked: that every panel falls inside the reported
    numerical range (widened by `tol`, in mm, for measurement scatter), and
    that the ordering the authors state -- single-step largest, detailed
    smallest -- comes out of the measurement.
    """
    lo, hi = PUBLISHED_RANGE["numerical"]
    rows = []
    for case in PUBLISHED:
        if keys and case["key"] not in keys:
            continue
        foil = load_foil(pdf, case["width"], case["height"], case["unit"],
                         pdf_page=case["page"], pdf_index=case["index"],
                         name=case["key"])
        res = measure(foil)
        if outdir:
            # show the authors' reported range on the histogram
            _write(res, outdir, case["key"], PUBLISHED_RANGE["numerical"])
        inside = lo - tol <= res.median <= hi + tol
        rows.append(dict(key=case["key"], rank=case["rank"],
                         measured=res.median, in_range=inside,
                         agreement=res.family_agreement(),
                         coherence=res.coherence, cite=case["cite"]))
        if verbose:
            print(f"  {case['key']:26s} {res.median:.2f} {case['unit']}  "
                  f"[IQR {res.stats()['q1']:.2f}-{res.stats()['q3']:.2f}]  "
                  f"families {res.family_agreement():3.0f} %  "
                  f"{'in range' if inside else 'OUTSIDE RANGE'}")
    ordered = [r["measured"] for r in sorted(rows, key=lambda r: r["rank"])]
    order_ok = all(a >= b - tol for a, b in zip(ordered, ordered[1:]))
    if verbose:
        print(f"  reported range {lo}-{hi} mm "
              f"({PUBLISHED_RANGE['source']}, p. {PUBLISHED_RANGE['page']})")
        print(f"  ordering single-step >= three-step >= detailed: "
              f"{'reproduced' if order_ok else 'NOT reproduced'}"
              + ("" if len(rows) == len(PUBLISHED) else "  (subset measured)"))
    passed = all(r["in_range"] for r in rows) and \
        (order_ok or len(rows) < len(PUBLISHED))
    return passed, rows
