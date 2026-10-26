# Hosting SootViz for free

The goal: a page a researcher can open, drop a soot foil onto, and get a cell
width from — at no running cost to us, and modelled on the Caltech Detonation
Database, which has stayed useful for twenty years because it is a static site
with no backend to rot.

## Recommendation: GitHub Pages + Pyodide

Run the Python in the visitor's browser. Pyodide is CPython compiled to
WebAssembly and ships NumPy, SciPy, Pillow and Matplotlib — which is the whole
of our dependency list. The site is then static files on GitHub Pages.

```
sootviz/                      (gh-pages)
  index.html                  drop target, controls, the verification step
  sootviz.js                  loads Pyodide, installs the wheel, calls measure()
  sootviz-0.1.0-py3-none-any.whl
  benchmarks/                 pre-rendered PNGs + JSON from sootviz.benchmark
  about.html                  method, credits, citation
```

Why this is the right choice:

- **Zero cost, zero maintenance.** No server, no container, no cold starts, no
  bill that arrives when a paper cites us and traffic spikes.
- **The images never leave the machine.** This matters more than it sounds:
  unpublished soot foils are data people are protective of, and "your image is
  not uploaded anywhere" removes the main reason a researcher would not try it.
  It also means no privacy policy, no data-retention question, no liability.
- **The same code runs in the browser and on the command line.** One package,
  one set of defaults, one validation suite. A browser-only reimplementation
  would drift from the Python and the published numbers would stop meaning
  anything.
- **It cannot go stale.** A static site keeps working after we stop touching it.

Costs to be honest about: Pyodide is a ~10 MB first load (cached afterwards),
and it runs perhaps 2-5x slower than native Python, so a 5 m plate takes a
minute or two rather than tens of seconds. Run the scale search in a web worker
so the page stays responsive, and show the per-stage progress the CLI prints.

## The verification step is part of the method, not decoration

The measurement can fail quietly: if the two families are misread, the number
is still produced and still looks plausible. So the page must not hand over a
number before the user has seen the measured crossings lying on the tracks.

Flow: drop image -> enter the physical extent -> measure -> **the verification
figure, full width, with "do the marks lie along the tracks?"** -> on "no",
the controls that matter (polarity, expected cell width, contrast window,
exclusion bands) with a short note on what each one does -> on "yes", the
number, the distribution, and the downloadable CSV/JSON.

The two internal checks go next to the number, in words: the two families are
independent measurements of the same quantity, so their disagreement is an
error bar, and the coherence says whether a two-family picture fits at all.

## Benchmarks on the site

Pre-rendered, never computed in the browser, so the page is evidence before the
visitor has uploaded anything:

1. **Synthetic foils with a known cell width** — the only true ground truth.
2. **Our KRAFT numerical foil** — the record from the paper.
3. **Published foils (Taileb et al. 2020)** — measured from the authors'
   own figures and compared with the range they report.

Each entry: the image, the verification figure, the distribution, the number,
and the reference value. Link the figures to the JSON so anyone can check.

Taileb et al.'s figures are copyrighted. The site shows *our* measurement
outputs and cites the paper; it does not redistribute their figures. Visitors
who want to reproduce that benchmark run it locally against their own copy of
the paper, which is what `sootviz benchmark --pdf` is for.

## Alternatives, and why not

| Option | Cost | Why not |
|---|---|---|
| Hugging Face Spaces (Gradio) | free CPU tier | images upload to a third party; sleeps; the free tier's terms can change |
| Streamlit Community Cloud | free | same upload concern, plus sleeping apps and a single-tenant resource limit |
| Google Colab notebook | free | great as a *secondary* entry point for people who want to script it; not a tool you can hand to a referee |
| A small VPS / Fly.io | ~$5/mo | a bill and a machine to patch, forever, for a tool that needs no server |

Worth doing anyway: ship a Colab notebook that `pip install`s the package, for
people who want to batch many foils, and publish the wheel to PyPI so
`pip install sootviz-compute` works. The browser page, the notebook and the CLI
then all run the same code.

## Citation and credit

Footer on every page, and in the README:

> SootViz — Karthik Thyagarajan, Strawhat CFD.
> Method described in Thyagarajan & Melguizo-Gavilanes (2026), Appendix C.

## Build steps (once the interface is wanted)

1. `python -m build` -> `sootviz-0.1.0-py3-none-any.whl`
2. `sootviz benchmark --out site/benchmarks` -> the gallery assets
3. `index.html` + `sootviz.js`: load Pyodide from the jsDelivr CDN,
   `micropip.install` the wheel, hand the dropped file to `load_foil`, render
   the figures Matplotlib produces into canvases.
4. Push to `gh-pages`; enable Pages on the repository.

Step 3 is the only new code, and it is a thin wrapper: everything it calls is
already in this package and already validated.
