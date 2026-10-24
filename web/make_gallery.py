"""Turn a `sootviz benchmark --out` directory into a small static gallery.

Run by CI after the benchmarks, so the site carries the figures behind the
numbers rather than only the numbers.
"""
import json
import sys
from pathlib import Path

HEAD = """<!doctype html>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>SootViz — benchmark gallery</title>
<link rel="stylesheet" href="../style.css">
<div class="wrap">
<header>
  <h1>Soot<span>Viz</span> — benchmarks</h1>
  <p class="tagline">The figures behind the numbers. Synthetic foils have a
  cell width known exactly by construction, so they are the only place where
  an error can be quoted against truth.</p>
  <p class="byline"><a href="../">back to the tool</a> ·
  <a href="https://strawhatcfd.com" target="_blank" rel="noopener">Strawhat CFD</a></p>
</header>
"""
TAIL = """<footer><p>Karthik Thyagarajan ·
<a href="https://strawhatcfd.com" target="_blank" rel="noopener">Strawhat CFD</a> ·
MIT licensed</p></footer>
</div>
"""


def main(d: Path) -> int:
    doc = json.loads((d / "benchmarks.json").read_text())
    out = [HEAD]
    for row in doc.get("synthetic", []):
        slug = "synthetic-" + row["name"].replace(", ", "-").replace(" ", "-")
        err = row["error_percent"]
        ok = "ok" if row["passed"] else "FAILED"
        out.append(f"""
<section class="step">
  <div class="step-head"><span class="num">{ok}</span><h2>{row['name']}</h2></div>
  <div class="stats">
    <div class="stat"><b>{row['truth']:.2f}</b><span>true λ (mm)</span></div>
    <div class="stat"><b>{row['measured']:.2f}</b><span>measured</span></div>
    <div class="stat"><b>{err:+.1f} %</b><span>error</span></div>
    <div class="stat"><b>{row['agreement']:.0f} %</b><span>families agree to</span></div>
    <div class="stat"><b>{row['coherence']:.2f}</b><span>coherence</span></div>
  </div>
  <figure><img src="{slug}_distribution.png" alt="distribution for {row['name']}">
    <figcaption>The result: every spacing measured against the true value, the
    two families separately, and the block medians along the record.</figcaption></figure>
  <figure><img src="{slug}_counting.png" alt="how spacings were counted for {row['name']}">
    <figcaption>Where the numbers come from — the vertical cuts, the
    trajectory crossings found on them, and the gaps between consecutive
    crossings, which are the measurements.</figcaption></figure>
  <figure><img src="{slug}_blocks.png" alt="blocks for {row['name']}">
    <figcaption>Each block with the family angles fitted inside it, its own
    cell width and its own count.</figcaption></figure>
  <figure><img src="{slug}_overview.png" alt="track families over {row['name']}">
    <figcaption>The per-block fields recombined over the whole record.</figcaption></figure>
  <figure><img src="{slug}_scale.png" alt="scale search for {row['name']}">
    <figcaption>The scale search. Flat means the answer does not depend on the
    assumed scale.</figcaption></figure>
</section>""")
    pub = doc.get("published", [])
    if pub:
        out.append('<h2>Published foils</h2><p style="font-size:14px">Measured '
                   'from the authors\' own figures. Those figures are '
                   'copyrighted and are not reproduced here; only our '
                   'measurements are shown.</p><table>'
                   '<tr><th>record</th><th class="num-col">SootViz</th>'
                   '<th>in reported range</th></tr>')
        for row in pub:
            out.append(f'<tr><td>{row["key"]}</td>'
                       f'<td class="num-col">{row["measured"]:.2f} mm</td>'
                       f'<td>{"yes" if row["in_range"] else "NO"}</td></tr>')
        out.append("</table>")
    out.append(TAIL)
    (d / "index.html").write_text("".join(out))
    print(f"gallery written: {d / 'index.html'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(Path(sys.argv[1] if len(sys.argv) > 1 else "_site/benchmarks")))
