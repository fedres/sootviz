"""Fast tests. The slow, physical validation lives in sootviz.benchmark."""
import numpy as np
import pytest

from sootviz import Config, measure, synth_foil, write_report
from sootviz.clean import clean, equalise
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
