// SootViz in the browser: Pyodide runs the same Python package the command
// line runs. Everything happens inside this worker, on the visitor's machine —
// no image is ever sent anywhere.
const PYODIDE = "https://cdn.jsdelivr.net/pyodide/v0.26.4/full/";
const WHEEL = "./sootviz_compute-0.1.0-py3-none-any.whl";

importScripts(PYODIDE + "pyodide.js");

let pyodide = null;

const say = (type, payload = {}) => self.postMessage({ type, ...payload });

async function boot() {
  say("stage", { text: "downloading the Python runtime (about 25 MB, cached afterwards)" });
  pyodide = await loadPyodide({ indexURL: PYODIDE });
  say("stage", { text: "loading numpy, scipy, matplotlib" });
  await pyodide.loadPackage(["numpy", "scipy", "matplotlib", "pillow", "micropip"]);
  say("stage", { text: "installing sootviz" });
  await pyodide.runPythonAsync(`
import micropip
await micropip.install("${WHEEL}")
import matplotlib
matplotlib.use("Agg")

import sys, js
from pyodide.ffi import to_js
class _Relay:
    """Send the pipeline's own progress lines to the page as they happen."""
    def write(self, t):
        t = t.rstrip()
        if t:
            js.postMessage(to_js({"type": "log", "text": t},
                                 dict_converter=js.Object.fromEntries))
    def flush(self):
        pass
sys.stdout = _Relay()
import sootviz
`);
  const v = pyodide.runPython("sootviz.__version__");
  say("ready", { version: v });
}

const bootPromise = boot().catch((e) => say("error", { text: String(e) }));

self.onmessage = async (ev) => {
  const m = ev.data;
  if (m.type !== "measure") return;
  await bootPromise;
  if (!pyodide) return;
  try {
    pyodide.FS.writeFile("/input_image", new Uint8Array(m.bytes));
    pyodide.globals.set("_opts", pyodide.toPy(m.opts));
    say("stage", { text: "measuring — the scale search evaluates the pipeline several times" });
    const res = await pyodide.runPythonAsync(`
import json, base64
from pathlib import Path
from sootviz import load_foil, measure, write_report
from sootviz.measure import Config

o = _opts.to_py() if hasattr(_opts, "to_py") else dict(_opts)
foil = load_foil("/input_image", float(o["width"]), float(o["height"]),
                 o.get("unit", "mm"), name=o.get("name", "foil"))
cfg = Config(
    cell=float(o["cell"]) if o.get("cell") else None,
    polarity=o.get("polarity", "auto"),
    grid=not o.get("no_grid", False),
    equalise_frac=float(o.get("equalise", 2.5)),
)
res = measure(foil, cfg, verbose=True)
written = write_report(res, Path("/report"), prefix="foil")

figs = {k: base64.b64encode(Path(v).read_bytes()).decode()
        for k, v in written.items() if str(v).endswith(".png")}
st = res.stats()
json.dumps(dict(
    figures=figs,
    median=st["median"], q1=st["q1"], q3=st["q3"], n=st["n"],
    mean=st["mean"], p5=st["p5"], p95=st["p95"],
    unit=foil.unit, cells_across=res.cells_across(),
    family_agreement=res.family_agreement(), coherence=res.coherence,
    sensitivity=res.scale_sensitivity(), polarity=res.polarity,
    blocks=len([b for b, _ in res.blocks if b is not None]),
    cell_scale=res.cell,
    report=Path("/report/foil_result.json").read_text(),
    spacings=Path("/report/foil_spacings.csv").read_text(),
))
`);
    say("result", { json: res });
  } catch (e) {
    say("error", { text: String(e) });
  }
};
