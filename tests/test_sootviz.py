"""Fast tests. The slow, physical validation lives in sootviz.benchmark."""
import numpy as np
import pytest

from sootviz import Config, measure, synth_foil, write_report
from sootviz.clean import clean, equalise
from sootviz.io import Foil
from sootviz.tracks import family_angles, grid_blocks, bank, ridge_map


# ---- the synthetic foil is the ground truth --------------------------------
def test_recovers_a_known_cell_width():
    f = synth_foil(width=24, height=16, s=30, lam=2.0, skew=0, jitter=0.05,
                   noise=0.05)
    r = measure(f, Config(cell=2.0 * 1.2, polarity="dark"))
    assert abs(r.median - 2.0) / 2.0 < 0.10
    assert r.family_agreement() < 15


def test_recovers_a_sheared_lattice():
    """A sheared pattern puts the two families at different angles; the cell
    width must still come out right, which is the point of separating them."""
    f = synth_foil(width=24, height=16, s=30, lam=2.0, skew=15, jitter=0.08)
    r = measure(f, Config(cell=2.4, polarity="dark"))
    assert abs(r.median - 2.0) / 2.0 < 0.15
    b = [x for x, _ in r.blocks if x is not None]
    assert b, "no block produced a pair of angles"
    assert max(abs(abs(x.a1) - abs(x.a2)) for x in b) > 5, \
        "the fitted families are symmetric, so the shear was not detected"


def test_scale_is_found_without_being_told():
    f = synth_foil(width=24, height=16, s=25, lam=2.0, skew=8, jitter=0.08)
    r = measure(f)                      # no cell given: bisection
    assert abs(r.median - 2.0) / 2.0 < 0.20


# ---- the pieces ------------------------------------------------------------
def test_grid_blocks_are_large_enough_to_carry_statistics():
    """A block must hold several cells, or its two modal angles are noise."""
    shape, s, cell = (400, 1200), 20.0, 1.0      # 20 x 60 units, 1-unit cells
    bl = grid_blocks(shape, s, cell)
    assert bl, "no blocks"
    w = bl[0][2]
    assert w >= 4 * cell - 1e-9
    assert w >= 0.5 * shape[0] / s - 1e-9
    assert len(bl) <= 60


def test_grid_blocks_never_explode_at_a_fine_trial_scale():
    bl = grid_blocks((400, 4000), 20.0, 0.02)
    assert len(bl) <= 60


def test_family_angles_rejects_noise():
    rng = np.random.default_rng(0)
    R = rng.random((160, 160))
    angles = np.arange(-60, 61, 10.0)
    _, _, score = family_angles(bank(R, 20.0, 1.0, angles), angles)
    assert score < 0.35, "pure noise was reported as a coherent two-family pattern"


def test_family_angles_finds_a_known_pair():
    f = synth_foil(width=20, height=20, s=25, lam=2.5, skew=0, jitter=0.0,
                   noise=0.0, vignette=0.0)
    R = ridge_map(f.image, f.s, 2.5)
    angles = np.arange(-70, 71, 5.0)
    a1, a2, score = family_angles(bank(R, f.s, 2.5, angles), angles)
    true = np.degrees(np.arctan(1 / 1.8))        # lam / (aspect * lam)
    assert min(abs(a1 - true), abs(a2 + true)) < 7.5
    assert score > 0.3


def test_equalise_preserves_the_pattern_it_is_given():
    f = synth_foil(width=30, height=20, s=25, lam=2.5, vignette=0.9)
    e = equalise(f.image, window_px=f.px(2.5 * 2.5))
    assert np.isfinite(e).all()
    # the vignette is a slow left-right trend: equalising must flatten it
    col = lambda a: a.mean(axis=0)
    span = lambda a: col(a).max() - col(a).min()
    assert span(e) / (e.std() + 1e-9) < span(f.image) / (f.image.std() + 1e-9)


def test_clean_runs_on_a_flat_image():
    out = clean(np.full((80, 80), 128.0), 20.0, 1.0)
    assert np.isfinite(out).all()


# ---- results and reporting -------------------------------------------------
def test_report_writes_what_it_promises(tmp_path):
    f = synth_foil(width=20, height=14, s=25, lam=2.0)
    r = measure(f, Config(cell=2.4, polarity="dark"))
    written = write_report(r, tmp_path, prefix="t")
    for role, p in written.items():
        assert p.exists() and p.stat().st_size > 0, f"{role} was not written"
    assert written["json"].suffix == ".json"
    assert written["csv"].suffix == ".csv"
    assert "verification" in written, \
        "the verification figure is the point of the report"


def test_statistics_are_consistent():
    f = synth_foil(width=20, height=14, s=25, lam=2.0)
    r = measure(f, Config(cell=2.4, polarity="dark"))
    st = r.stats()
    assert st["q1"] <= st["median"] <= st["q3"]
    assert st["n"] == len(r.lam)
    assert r.cells_across() == pytest.approx(f.height / r.median)


def test_empty_image_fails_loudly():
    f = synth_foil(width=10, height=10, s=20, lam=2.0, contrast=0.0,
                   noise=0.0, vignette=0.0)
    f.image[:] = 128.0
    with pytest.raises(RuntimeError):
        measure(f, Config(polarity="dark"))


# ---- the block-level work, and the evidence behind the count ---------------
def test_every_spacing_comes_from_a_recorded_cut():
    """The counting figure is only honest if the cuts it draws are the cuts the
    statistics were taken from, so the result must carry them."""
    f = synth_foil(width=30, height=20, s=25, lam=2.4, skew=10, jitter=0.1)
    r = measure(f, Config(cell=2.9, polarity="dark"))
    assert r.cuts1 and r.cuts2, "no cuts recorded"
    for x, peaks in r.cuts1[:5]:
        assert 0 <= x <= f.width
        assert np.all((peaks >= 0) & (peaks <= f.height))
    # a cut with k peaks contributes k-1 spacings, so the totals must agree
    expect = sum(max(0, len(p) - 1) for _, p in r.cuts1)
    assert len(r.lam1) == expect, "the histogram does not match the drawn cuts"
    assert len(r.lam_x) == len(r.lam)


def test_blocks_carry_their_own_angles_and_measurement():
    f = synth_foil(width=40, height=20, s=25, lam=2.4, skew=12, jitter=0.1)
    r = measure(f, Config(cell=2.9, polarity="dark"))
    rows = r.per_block()
    assert len(rows) == len([b for b, _ in r.blocks if b is not None])
    assert len(rows) > 1, "a 40x20 record should be split into several blocks"
    for b in rows:
        assert b["a1"] > 0 > b["a2"], "the two families must straddle zero"
        assert b["n"] > 0 and np.isfinite(b["lam"])
    # blocks tile the record left to right
    xs = [b["x"] for b in rows]
    assert min(xs) < f.width / 3 < max(xs)


def test_report_writes_the_diagnostics(tmp_path):
    f = synth_foil(width=30, height=20, s=25, lam=2.4, skew=10)
    r = measure(f, Config(cell=2.9, polarity="dark"))
    w = write_report(r, tmp_path, prefix="d")
    for key in ("distribution", "counting", "blocks", "overview", "clean",
                "verification", "csv", "blocks_csv", "json"):
        assert key in w, f"{key} missing from the report"
        assert w[key].exists() and w[key].stat().st_size > 0
    # the scale panel plots the search, so it exists only when one was run
    assert ("scale" in w) == bool(r.ladder)
    head, *body = w["blocks_csv"].read_text().strip().splitlines()
    assert head.startswith("x[") and "angle1[deg]" in head
    assert len(body) == len(r.per_block())
    spac = w["csv"].read_text().strip().splitlines()
    assert spac[0].endswith("family")
    assert len(spac) - 1 == len(r.lam)


def test_result_json_is_strictly_valid(tmp_path):
    """NaN is not JSON. A run given its scale has no scale search, so its
    sensitivity is undefined -- and Python would write a bare NaN, which every
    strict parser rejects, discarding the whole result rather than one field."""
    import json
    f = synth_foil(width=24, height=16, s=25, lam=2.4, skew=10)
    r = measure(f, Config(cell=2.9, polarity="dark"))
    assert not np.isfinite(r.scale_sensitivity()), "this test needs the NaN case"
    text = write_report(r, tmp_path, prefix="j")["json"].read_text()
    for token in ("NaN", "Infinity", "-Infinity"):
        assert token not in text, f"{token} would break a strict JSON parser"
    doc = json.loads(text)                      # strict by default in Python
    assert doc["scale_sensitivity_percent"] is None
    assert doc["stats"]["median"] > 0


# ---- is the measurement following the record, or its own filters? ----------
def test_crossings_land_on_ridges_not_at_chance():
    """The alignment check is measured against the ridge map, which the cuts
    never see -- they are taken on the oriented family fields -- so it is an
    independent check that the measurement is following the record."""
    f = synth_foil(width=30, height=20, s=25, lam=2.5, skew=10, jitter=0.1,
                   noise=0.25)
    r = measure(f)
    assert r.alignment() > 75, f"crossings near chance: {r.alignment():.0f}"
    assert r.on_track() > 40
    assert len(r.align) == sum(len(p) for _, p in r.cuts1) \
        + sum(len(p) for _, p in r.cuts2)


def test_noise_alone_scores_near_chance():
    """A record with no cellular structure must not produce a confident
    number: the crossings should sit near the 50th percentile."""
    rng = np.random.default_rng(3)
    f = synth_foil(width=20, height=14, s=25, lam=2.0, contrast=0.0,
                   noise=0.0, vignette=0.0)
    f.image[:] = 128 + 30 * rng.normal(size=f.image.shape)
    try:
        r = measure(f, Config(cell=2.4, polarity="dark"))
    except RuntimeError:
        return                      # refusing outright is also a pass
    assert r.alignment() < 80 or r.coherence < 0.3, (
        f"pure noise scored alignment {r.alignment():.0f}, coherence "
        f"{r.coherence:.2f} -- the checks are not discriminating")


def test_answer_does_not_depend_on_scan_resolution():
    """Every filter is a fraction of a cell, so the analysis is done at a fixed
    cells-per-pixel. The same record scanned coarsely or finely must agree."""
    a = synth_foil(width=30, height=20, s=18, lam=2.5, skew=10, jitter=0.08)
    b = synth_foil(width=30, height=20, s=54, lam=2.5, skew=10, jitter=0.08)
    ra, rb = measure(a), measure(b)
    assert abs(ra.median - rb.median) / rb.median < 0.10, (
        f"{ra.median:.3f} at {a.s:.0f} px/mm vs {rb.median:.3f} at "
        f"{b.s:.0f} px/mm")


def test_a_large_record_does_not_exhaust_memory():
    """The filter bank holds one response per angle over the whole frame, so
    the frame has to be bounded or a big photograph fails with a memory error
    rather than a message."""
    from sootviz.measure import _decimate
    f = synth_foil(width=60, height=12, s=120, lam=2.0)      # 8.6 Mpx
    assert f.image.size > 8e6
    d = _decimate(f, 2.0, Config().px_per_cell)
    assert d.image.size <= 1.6e6, f"frame not capped: {d.image.size/1e6:.1f} Mpx"


# ---- claims the paper makes about determinism ------------------------------
def test_measurement_is_deterministic():
    """The paper states the pipeline is a deterministic function of the record
    and its extent. Same input, same output, bit for bit."""
    f = synth_foil(width=24, height=16, s=25, lam=2.4, skew=10, jitter=0.1)
    a = measure(f, Config(cell=2.9, polarity="dark"))
    b = measure(f, Config(cell=2.9, polarity="dark"))
    assert a.median == b.median
    assert np.array_equal(np.sort(a.lam), np.sort(b.lam))
    assert a.cell == b.cell and a.coherence == b.coherence


def test_answer_is_independent_of_the_unit_used():
    """Every length in the pipeline is a fraction of the cell width, so the
    same record declared in different units must give the same physical
    answer."""
    img = synth_foil(width=30, height=20, s=25, lam=2.5, skew=8,
                     jitter=0.08).image
    mm = Foil(img, 25.0, 30.0, 20.0, "mm", "mm")
    cm = Foil(img, 250.0, 3.0, 2.0, "cm", "cm")   # same record, declared in cm
    a = measure(mm, Config(polarity="dark"))
    b = measure(cm, Config(polarity="dark"))
    assert abs(a.median - 10 * b.median) / a.median < 0.02, (
        f"{a.median:.4f} mm vs {b.median:.4f} cm")
