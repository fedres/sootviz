# SootViz-compute

Measuring the detonation cell width from a soot-foil record — numerical or
photographed — with the measurement you can check by eye before you trust it.

*Karthik Thyagarajan · [Strawhat CFD](https://strawhatcfd.com)*

**[Open the web version →](https://fedres.github.io/sootviz/)** — runs in your browser, your image never leaves your machine.

---

## Why another cell-size tool

A soot foil carries two families of triple-point trajectories, and the cell
width is the **transverse spacing between adjacent trajectories of one
family**. Two things make that hard to get right, and both are easy to get
wrong silently:

- **The pattern is usually sheared.** The two families sit at different angles,
  so the trajectory bounding a cell transversely is *not* the one directly
  above it. A nearest-site-above rule produces a long tail of spurious wide
  cells. (This is a mistake we made and had to go back and fix.)
- **Autocorrelation and spectra measure the wrong thing.** On a real foil the
  cells group into larger patches, and the dominant periodicity of the raw
  image is the *grouping*, not the cell. On one experimental plate that error
  was a factor of three.

So SootViz separates the two families first, measures each one separately, and
uses their agreement as an internal check.

## The method

1. **Clean** — despeckle, then equalise the contrast over a window several
   cells wide. The window matters: equalising at the scale of a cell erases the
   pattern you are trying to measure.
2. **Ridge map** — a difference of Gaussians, signed by whether tracks are dark
   (numerical foils, usually max-pressure plots) or bright (scoured soot).
3. **Oriented filter bank** — line matched filters from −70° to +70°. Each
   smooths along its own direction and not across it, so a trajectory responds
   only to the filters aligned with it. The two modal angles are the two
   families.
4. **Local angles, global measurement** — on a long record the pattern turns,
   so the angles are fitted per block and the per-block fields are blended with
   a partition of unity into two continuous fields. The spacings are then read
   off full-height cuts through those fields: measuring inside a block biases
   the result low, because peak prominence is judged against each cut's own
   maximum and a short cut makes weak maxima look prominent.
5. **Spacings** — peaks along each vertical cut; consecutive differences are
   cell widths. Each family gives an independent estimate.

Every length in the pipeline is a fraction of the cell width — which is the
quantity being measured, so the scale has to be chosen rather than assumed.
Measured against the assumed scale c, the result λ̂(c) has two regimes: once the
filters are fine enough to separate adjacent trajectories, **λ̂ is flat** — the
answer stops depending on the assumption, which is the signature of being
right — and below that, filters finer than the pattern begin subdividing real
spacings and λ̂ follows c down. Coming from coarse, SootViz takes the knee: the
coarsest scale whose measured spacing is consistent with the scale used to
measure it. Every run writes this curve out (`*_scale.png`), so you can see
whether your measurement sits on a plateau or on a slope.

Two simpler rules were tried first and are both wrong, which is worth
recording. Fixed-point iteration stops wherever its own tolerance happens to
bite, and on a noisy synthetic foil it collapsed to 1.8 mm for a true 3.0 mm.
Root-finding on λ̂(c) − c looks more principled, but on a regular lattice there
is no root near the answer — λ̂ stays flat at 2.5 mm while c falls past it — and
the only root lies deep in the subdividing regime, returning 1.5 mm for a
lattice built with 2.5 mm cells.

## Two checks that come free

- **Family agreement.** The two families are independent measurements of the
  same quantity. Their difference is an error bar, not a diagnostic to be
  dismissed — above about 15 % the geometry has not been read correctly.
- **Coherence.** The fraction of strong track pixels belonging to the two modal
  directions. Below about 0.25 a two-family picture does not fit the image at
  all, and no number should be quoted.

Both are printed with every result and written into the JSON.

## What the number means

For a regular pattern λ̂(c) is genuinely flat and the cell width is a sharp
quantity. For an irregular one it is not: an irregular cellular front carries
structure over a range of scales, the measured λ depends on the scale you look
at, and that is physics rather than a defect of the method — it is why
published cell widths for irregular mixtures differ between authors by factors
approaching two. So:

- quote the **distribution**, not only the median — the IQR is written into
  every report, and `*_spacings.csv` holds every spacing measured;
- read the **scale panel**. A flat curve means the number is robust; a sloping
  one means the record has no single cell width, and the number you quote
  should be accompanied by the scale at which it was taken;
- compare like with like. Two foils measured by this pipeline with the same
  settings can be compared even where the absolute value is scale-dependent,
  because the bias is common to both.

## Install

```bash
pip install -e .                 # or: pip install -r requirements.txt
pip install -e '.[pdf]'          # to read foils straight out of a PDF figure
```

## Use

```bash
# a photographed or exported foil: the physical extent is required
sootviz measure foil.png --width 120 --height 38 --unit mm --out run1

# a figure inside a PDF, third raster on page 19
sootviz measure paper.pdf --width 30 --height 20 --pdf-page 19 --pdf-index 2

# what rasters does this PDF page hold?
sootviz rasters paper.pdf --pdf-page 19

# a synthetic foil with a known cell width, to see the pipeline work
sootviz synth --lam 2.5 --skew 12 --noise 0.3 --out demo

# the validation suite
sootviz benchmark                      # synthetic only
sootviz benchmark --pdf CNFTaileb.pdf  # plus published foils
```

Then **look at `*_families.png` before you use the number.** The contours must
lie along the tracks, one colour per family. If they do not, the controls that
matter are, in order: `--polarity`, `--cell` (an expected cell width, which
turns the scale search off), `--equalise` (the contrast window, in cells), and
`--exclude-x` (physical x positions to skip — seams, mounting damage, burn-through).

From Python, see [examples/python_api.py](examples/python_api.py).

## Validation

Four levels of evidence, in increasing order of realism. Reproduce with
`sootviz benchmark` (the published level needs your own copy of the paper).

**1. Synthetic foils — the only place truth is known exactly.**

| case | true λ | measured | error | families | coherence |
|---|---|---|---|---|---|
| regular | 2.50 mm | 2.48 | −0.8 % | 0 % | 0.84 |
| sheared (asymmetric families) | 2.50 mm | 2.48 | −0.8 % | 2 % | 0.71 |
| irregular (0.20 λ jitter) | 2.00 mm | 1.88 | −6.0 % | 0 % | 0.56 |
| noisy, low contrast | 3.00 mm | 2.52 | −16.0 % | 6 % | 0.50 |

The last case is the documented limit: heavy noise with 0.35 contrast shortens
the plateau, and the scale search — which deliberately prefers the finer end,
because on real records the danger is drifting from the cell onto the
*grouping* of cells — then lands about 15 % low. A real record in that state
reports a high scale sensitivity, which is the signal to quote a distribution
rather than a number.

**2. Our own numerical foil** (KRAFT, one-step chemistry, 30 × 20 mm):
λ = **2.35 mm**, IQR 1.79–3.06, n = 1345, families agree to 10 %, scale
sensitivity 6 %. An earlier hand measurement of the same record with a single
absolute filter gave 2.28 mm (+3.1 %).

**3. Published foils — Taileb, Melguizo-Gavilanes & Chinnayya,**
*Combust. Flame* **218** (2020) 247, Fig. 5. The authors report a range over
their three chemistry models rather than a value per panel — "the cell size was
numerically estimated to lie in the range of 0.9 mm ≤ λ_num ≤ 1.9 mm" (p. 17)
— and state that detailed chemistry gives the smaller cell, while single-step
was tuned to the experimental size. So what can honestly be checked is range
containment and that ordering:

| panel (their chemistry) | SootViz | IQR | in reported range | 
|---|---|---|---|
| single-step | 2.07 mm | 1.50–2.74 | yes |
| three-step chain-branching | 1.28 mm | 0.93–1.72 | yes |
| detailed | 0.93 mm | 0.71–1.19 | yes |

Ordering single-step ≥ three-step ≥ detailed: **reproduced**, and the two
extremes land on the authors' two bounds. Their figures are copyrighted and are
not redistributed here; point `--pdf` at your own copy.

**4. A photographed experimental plate**, 5 m × 1 m, where the pattern turns
along the record: λ = **0.162 m**, 13 blocks, family angles swinging from
+35°…+65° to −65°…−45° along the plate, per-block λ falling 0.18 → 0.14 m
downstream. This record is the instructive one: the families disagree by 19 %
and the scale sensitivity is 23 %, so **SootViz warns on both counts** rather
than handing over a clean-looking number. The honest reading is λ ≈ 0.16–0.24 m
depending on the scale examined. For contrast, autocorrelation of the same
image returns 0.52–0.59 m — it locks onto the grouping of cells, not the cell,
and is wrong by a factor of three.


## Hosting this as a web tool

See [docs/HOSTING.md](docs/HOSTING.md). The short version: GitHub Pages plus
Pyodide, so the same Python runs in the visitor's browser, **the image never
leaves their machine**, and there is no server to pay for or maintain — a static
site in the spirit of the Caltech Detonation Database.

## Using it on the web

A browser version runs the same Python through WebAssembly, so **the image
never leaves the machine it is opened on** — nothing is uploaded, logged or
stored:

**https://fedres.github.io/sootviz/**

It walks through loading a foil, giving its physical size, measuring, and —
before any number is shown — the verification overlay and the question *do the
contours follow the tracks?* The benchmark gallery behind the numbers below is
at `/benchmarks/`.

## If you use SootViz, please cite it

It is free and always will be. Citing it is what makes it worth maintaining,
and it lets a reader check how your cell widths were obtained.

```bibtex
@software{sootviz,
  author  = {Thyagarajan, Karthik},
  title   = {{SootViz}: cell-width measurement from detonation soot-foil records},
  year    = {2026},
  note    = {Strawhat CFD},
  url     = {https://github.com/fedres/sootviz}
}
```

and the paper describing the method (Appendix C):

```bibtex
@article{thyagarajan2026lbfs,
  author  = {Thyagarajan, Karthik and Melguizo-Gavilanes, Jos\'e},
  title   = {A systematic assessment of the lattice {B}oltzmann flux solver with
             face-local gating as a scheme for flame acceleration, detonation
             and deflagration-to-detonation transition. {P}art {I}:
             one-dimensional detonation dynamics and inviscid detonations},
  year    = {2026}
}
```

In text: *cell widths were measured with SootViz (Thyagarajan, Strawhat CFD),
which separates the two triple-point trajectory families and reports their
disagreement as an error bar.*

## Who made this

Built by **Karthik Thyagarajan** at **[Strawhat CFD](https://strawhatcfd.com)**.
Questions, bug reports and foils that defeat it are all welcome — open an issue.

## Licence

MIT — see [LICENSE](LICENSE).
