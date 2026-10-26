# SootViz-compute

Measuring the detonation cell width from a soot-foil record — numerical or
photographed — with the measurement you can check by eye before you trust it.

*Karthik Thyagarajan · [Strawhat CFD](https://strawhatcfd.com)*

**[Open the web version →](https://fedres.github.io/sootviz/)** — runs in your
browser, your image never leaves your machine.
**[See the benchmark gallery →](https://fedres.github.io/sootviz/benchmarks/)** —
the figures behind every number below.

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
   pattern you are trying to measure. It is chosen from the record, by trying
   several and keeping the one whose two families come out most coherent.
2. **Ridge map** — a difference of Gaussians, signed by whether tracks are dark
   (numerical foils, usually max-pressure plots) or bright (scoured soot).
   Which one is read from the record: trajectories are thin and what lies
   between them is the wide interior of a cell, so the correct sign gives a
   concentrated, heavy-tailed ridge map and the wrong one spreads its response
   over the interiors. Nothing needs to be set by hand.
3. **Oriented filter bank** — line matched filters from −70° to +70°, with
   angles in image coordinates, so a positive angle is a line that descends to
   the right as you look at it. Each smooths along its own direction and not
   across it, so a trajectory responds only to the filters aligned with it.
   The two modal angles are the two families.
4. **Local angles, global measurement** — on a long record the pattern turns,
   so the angles are fitted per block and the per-block fields are blended with
   a partition of unity into two continuous fields. The spacings are then read
   off full-height cuts through those fields: measuring inside a block biases
   the result low, because peak prominence is judged against each cut's own
   maximum and a short cut makes weak maxima look prominent.
5. **Spacings** — peaks along each vertical cut; consecutive differences are
   cell widths. Each family gives an independent estimate.

The analysis always happens at a fixed number of pixels per cell, so the same
record measured from a small scan and a large one gives the same answer — and
the cross-track filter, the thing that selects a thin track, is never allowed
to fall below a pixel and stop filtering.

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
- **Alignment.** For every crossing the measurement found, the percentile rank
  of ridge strength at that point within its own column. 50 is chance; around
  90 is a clean record, and the high seventies still means the cuts are
  finding soot tracks rather than their own filters. This is checked against
  the ridge map, which the cuts never see — they are taken on the oriented
  family fields — so it is an independent check and not a restatement of how
  the peaks were chosen.

Both are printed with every result and written into the JSON.

## What the number means

For a regular pattern λ̂(c) is genuinely flat and the cell width is a sharp
quantity. For an irregular one it is not: an irregular cellular front carries
structure over a range of scales, the measured λ depends on the scale you look
at, and that is physics rather than a defect of the method — it is why
published cell widths for irregular mixtures differ between authors by factors
approaching two. So:

- quote the **distribution**, not only the median — `*_distribution.png` shows
  it three ways and `*_spacings.csv` holds every spacing measured, tagged with
  its position and family;
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
# a photographed or exported foil: the physical extent is all that is required.
# Whether the tracks are dark or bright, the contrast window and the cell scale
# are all read from the record.
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

### What a run writes

| file | what it is |
|---|---|
| `*_distribution.png` | **the result**: every spacing, the two families separately, and the block medians along the record |
| `*_families.png` | the verification overlay, in strips — check this first |
| `*_counting.png` | the cuts and the detected trajectory crossings, with spacings labelled. Every number in the histogram is one gap here |
| `*_blocks.png` | each block with its own fitted family angles, cell width and count |
| `*_overview.png` | the per-block fields recombined over the whole record |
| `*_clean.png` | before and after despeckling and local contrast equalisation |
| `*_scale.png` | the scale search, so you can see whether you are on a plateau |
| `*_spacings.csv` | every spacing, with its x position and which family it came from |
| `*_blocks.csv` | per-block angles, cell width, count and coherence |
| `*_result.json` | everything above as data, including the scale search |

Then **look at `*_families.png` before you use the number.** Each mark is one
trajectory crossing one cut, and every spacing in the histogram is the gap
between two marks of the same colour on the same cut — so the marks should lie
along the tracks and string out into lines following them. If they do not, the
controls that matter are, in order: `--polarity` (normally read from the
record), `--cell` (an expected cell width, which turns the scale search off),
`--equalise` (the contrast window in cells, normally chosen by trying several
and keeping the most coherent), and `--exclude-x` (physical x positions to
skip — seams, mounting damage, burn-through).

From Python, see [examples/python_api.py](examples/python_api.py).

## Validation

Four levels of evidence, in increasing order of realism. Reproduce with
`sootviz benchmark` (the published level needs your own copy of the paper).

**1. Synthetic foils — the only place truth is known exactly.**

| case | true λ | measured | error | families | alignment |
|---|---|---|---|---|---|
| regular | 2.50 mm | 2.49 | −0.3 % | 1 % | 93 |
| sheared (asymmetric families) | 2.50 mm | 2.49 | −0.3 % | 1 % | 92 |
| irregular (0.20 λ jitter) | 2.00 mm | 1.93 | −3.7 % | 0 % | 89 |
| noisy, low contrast | 3.00 mm | 2.55 | −15.0 % | 2 % | 90 |

**2. The same lattice, degraded the way real plates are.** A known 3.00 mm
lattice, each defect applied on its own:

| degradation | measured | error | alignment | families |
|---|---|---|---|---|
| clean | 3.01 mm | +0.4 % | 93 | 0 % |
| uneven lighting (strong vignette) | 3.01 mm | +0.4 % | 92 | 1 % |
| heavy soot grain | 2.98 mm | −0.8 % | 91 | 0 % |
| faint tracks (0.22 contrast) | 3.01 mm | +0.4 % | 93 | 0 % |
| irregular front (0.30 λ jitter) | 2.83 mm | −5.5 % | 89 | 1 % |
| all of them at once | 2.52 mm | −16.0 % | 90 | 12 % |

Lighting, grain and faint contrast cost essentially nothing on their own —
that is what the cleaning stage is for. Only genuine irregularity moves the
answer, and only the pile-up of every defect at once costs more than a few per
cent. Note that **alignment stays near 90 throughout**: even where the answer
drifts, the measurement is still following the record rather than inventing
structure, which is exactly what that check is for.

**3. Our own numerical foil** (KRAFT, one-step chemistry, 30 × 20 mm):
λ = **2.39 mm**, IQR 1.79–3.10, n = 1233, families agree to 11 %, alignment 90
with 83 % of crossings in the top quartile, scale sensitivity 4 %.

**4. Published foils — Taileb, Melguizo-Gavilanes & Chinnayya,**
*Combust. Flame* **218** (2020) 247, Fig. 5. The authors report a range over
their three chemistry models rather than a value per panel — "the cell size was
numerically estimated to lie in the range of 0.9 mm ≤ λ_num ≤ 1.9 mm" (p. 17)
— and state that detailed chemistry gives the smaller cell, while single-step
was tuned to the experimental size. So what can honestly be checked is range
containment and that ordering:

| panel (their chemistry) | SootViz | IQR | in reported range |
|---|---|---|---|
| single-step | 2.02 mm | 1.49–2.64 | yes |
| three-step chain-branching | 1.22 mm | 0.88–1.66 | yes |
| detailed | 0.88 mm | 0.66–1.16 | yes |

Ordering single-step ≥ three-step ≥ detailed: **reproduced**, and the two
extremes land on the authors' two bounds. Their figures are copyrighted and are
not redistributed here; point `--pdf` at your own copy.

**5. A photographed experimental plate**, 5 m × 1 m, where the pattern turns
along the record: λ = **0.181 m** over 13 blocks, family angles swinging from
+35°…+65° to −65°…−45°, per-block λ falling downstream. The two families agree
to 2 % and the crossings sit at the 83rd percentile of ridge strength. Scale
sensitivity is 23 %, so this record genuinely has no single cell width and the
distribution is the answer — the histogram is visibly **bimodal**, with
populations near 0.12 and 0.19 m that the median falls between. For contrast,
autocorrelation of the same image returns 0.52–0.59 m: it locks onto the
grouping of cells, not the cell, and is wrong by a factor of three.

Measuring the same plate from a 16× larger scan gives 0.184 m — the analysis
always happens at the same cells-per-pixel, so the answer does not depend on
how finely the record was scanned.

## Using it on the web
## Using it on the web

A browser version runs the same Python through WebAssembly, so **the image
never leaves the machine it is opened on** — nothing is uploaded, logged or
stored:

**https://fedres.github.io/sootviz/**

It walks through loading a foil, giving its physical size, measuring, and —
before any number is shown — the verification overlay and the question *do the
marks lie along the tracks?* The benchmark gallery behind the numbers below is
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
