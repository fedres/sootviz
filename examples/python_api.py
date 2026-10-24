"""Measuring one foil from Python, and reading the result.

    python examples/python_api.py
"""
from pathlib import Path

from sootviz import Config, load_foil, measure, synth_foil, write_report

# A real foil: give the physical extent, because every length in the pipeline
# is physical. An image alone cannot be measured.
#   foil = load_foil("myfoil.png", width=120.0, height=38.0, unit="mm")
# Straight out of a PDF figure (needs the `pdf` extra):
#   foil = load_foil("paper.pdf", 30.0, 20.0, pdf_page=19, pdf_index=0)

foil = synth_foil(width=40, height=25, s=25, lam=3.0, skew=10, jitter=0.12)

res = measure(foil, Config(), verbose=True)
print(res.summary())

# LOOK AT THIS FIGURE before using the number: the contours must lie along the
# tracks, one colour per family. Everything else follows from that being right.
out = write_report(res, Path("sootviz_out"), reference=3.0, prefix="example")
for role, path in out.items():
    print(f"  {role:13s} {path}")

st = res.stats()
print(f"\nlambda = {st['median']:.3g} {foil.unit}  "
      f"(IQR {st['q1']:.3g}-{st['q3']:.3g}, n = {st['n']})")
print(f"the two families differ by {res.family_agreement():.0f} % "
      f"-- treat that as the error bar")
