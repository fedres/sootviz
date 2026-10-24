// SootViz page logic. The measurement itself lives in worker.js, which runs
// the Python package through Pyodide; this file only drives the interface.
const $ = (id) => document.getElementById(id);
const show = (id, on = true) => { $(id).hidden = !on; };

let bytes = null, result = null, worker = null;

// ── loading a foil ──────────────────────────────────────────────────────────
const drop = $("drop"), file = $("file");
drop.addEventListener("click", () => file.click());
drop.addEventListener("dragover", (e) => { e.preventDefault(); drop.classList.add("over"); });
drop.addEventListener("dragleave", () => drop.classList.remove("over"));
drop.addEventListener("drop", (e) => {
  e.preventDefault(); drop.classList.remove("over");
  if (e.dataTransfer.files[0]) take(e.dataTransfer.files[0]);
});
file.addEventListener("change", () => file.files[0] && take(file.files[0]));

async function take(f) {
  if (!f.type.startsWith("image/")) {
    hint("That is not an image. For a figure inside a PDF, use the command line tool.");
    return;
  }
  bytes = await f.arrayBuffer();
  $("preview-img").src = URL.createObjectURL(f);
  show("preview");
  drop.innerHTML = `<strong>${f.name}</strong><small>${(f.size / 1024).toFixed(0)} kB — click to choose a different one</small>`;
  ready();
}

const hint = (t) => { $("hint").textContent = t; };
function ready() {
  const ok = bytes && parseFloat($("w").value) > 0 && parseFloat($("h").value) > 0;
  $("run").disabled = !ok;
  hint(ok ? "" : bytes ? "Give the physical width and height of the image."
                       : "Load an image and give its size.");
}
["w", "h"].forEach((id) => $(id).addEventListener("input", ready));

// ── running ─────────────────────────────────────────────────────────────────
$("run").addEventListener("click", () => {
  show("step-run"); show("step-check", false); show("step-result", false);
  $("log").textContent = ""; $("run").disabled = true;
  $("stage").textContent = "starting…";
  if (!worker) {
    worker = new Worker("worker.js");
    worker.onmessage = handle;
    worker.onerror = (e) => fail(e.message || "the worker failed to start");
  }
  worker.postMessage({
    type: "measure",
    bytes,
    opts: {
      width: parseFloat($("w").value), height: parseFloat($("h").value),
      unit: $("u").value, cell: parseFloat($("cell").value) || 0,
      polarity: $("pol").value, equalise: parseFloat($("eq").value) || 2.5,
      name: "your foil",
    },
    // structured-cloned, not transferred: the buffer stays usable here so the
    // same foil can be measured again after adjusting the controls.
  });
});

function handle(ev) {
  const m = ev.data;
  if (m.type === "stage") { $("stage").textContent = m.text; }
  else if (m.type === "log") {
    const el = $("log");
    el.textContent += m.text + "\n";
    el.scrollTop = el.scrollHeight;
  } else if (m.type === "result") { arrived(JSON.parse(m.json)); }
  else if (m.type === "error") { fail(m.text); }
}

function fail(text) {
  $("stage").textContent = "It did not work.";
  $("log").textContent += "\n" + text + "\n";
  $("run").disabled = false;
}

const img = (id, b64) => { if (b64) $(id).src = "data:image/png;base64," + b64; };

function arrived(r) {
  result = r;
  $("stage").textContent = "done — now check it.";
  img("fig-verify", r.figures.overview || r.figures.verification);
  img("fig-scale", r.figures.scale);
  img("fig-dist", r.figures.distribution);
  show("step-check");
  $("run").disabled = false;
  $("step-check").scrollIntoView({ behavior: "smooth", block: "start" });
}

// ── the approval gate ───────────────────────────────────────────────────────
$("no").addEventListener("click", () => {
  show("step-check", false); show("step-result", false);
  $("step-load").scrollIntoView({ behavior: "smooth", block: "start" });
  hint("Try, in this order: set how the tracks appear; give an expected λ near "
     + "the size you can see; raise the contrast window. Then measure again.");
  $("pol").focus();
});

$("yes").addEventListener("click", () => {
  const r = result, u = r.unit;
  const sig = (x, n = 3) => Number(x).toPrecision(n).replace(/\.?0+$/, "");
  $("lambda").innerHTML = `λ = ${sig(r.median)} ${u} <small>— the median spacing</small>`;
  $("stats").innerHTML = [
    ["interquartile range", `${sig(r.q1)}–${sig(r.q3)} ${u}`],
    ["spacings measured", r.n],
    ["cells across the height", sig(r.cells_across, 2)],
    ["families agree to", `${Math.round(r.family_agreement)} %`],
    ["scale sensitivity", isFinite(r.sensitivity) ? `${Math.round(r.sensitivity)} %` : "—"],
    ["coherence", sig(r.coherence, 2)],
    ["tracks read as", r.polarity],
  ].map(([k, v]) => `<div class="stat"><b>${v}</b><span>${k}</span></div>`).join("");

  const flags = [];
  if (r.family_agreement > 15)
    flags.push([`The two families disagree by ${Math.round(r.family_agreement)} %.`
      + " They are independent measurements of the same quantity, so this is a"
      + " real error bar on the number — above about 15 % the geometry may have"
      + " been read wrongly in part of the record.", "bad"]);
  if (r.coherence < 0.25)
    flags.push([`Orientation coherence is ${sig(r.coherence, 2)}. Below about`
      + " 0.25 a two-family picture does not fit this image, and no single cell"
      + " width should be quoted from it.", "bad"]);
  if (isFinite(r.sensitivity) && r.sensitivity > 25)
    flags.push([`The answer moves ${Math.round(r.sensitivity)} % across`
      + " neighbouring scales, so this record has no single cell width. That is"
      + " physics, not a fault: quote the distribution and the scale together.", ""]);
  if (!flags.length)
    flags.push(["Both internal checks are clean: the families agree and the"
      + " answer sits on a plateau in scale.", "ok"]);
  $("flags").innerHTML = flags.map(([t, k]) =>
    `<div class="flag ${k === "ok" ? "" : k}" ${k === "ok" ? 'style="border-left-color:var(--ok)"' : ""}>${t}</div>`).join("");

  const blob = (text, type) => URL.createObjectURL(new Blob([text], { type }));
  $("dl-csv").href = blob(r.spacings, "text/csv");
  $("dl-json").href = blob(r.report, "application/json");

  show("step-result");
  $("step-result").scrollIntoView({ behavior: "smooth", block: "start" });
});

$("again").addEventListener("click", () => {
  show("step-result", false); show("step-check", false); show("step-run", false);
  $("step-load").scrollIntoView({ behavior: "smooth", block: "start" });
});

// ── citation ────────────────────────────────────────────────────────────────
$("copy-cite").addEventListener("click", async (e) => {
  await navigator.clipboard.writeText($("bib").textContent + "\n\n" + $("bib2").textContent);
  e.target.textContent = "Copied";
  setTimeout(() => (e.target.textContent = "Copy both"), 1600);
});
