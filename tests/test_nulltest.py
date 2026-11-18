import numpy as np
import pytest

from polynome2pi.evaluation.nulltest import (
    Target,
    correction_series,
    evaluate_target,
    local_density,
    tuning_capacity,
)


@pytest.fixture
def uniform_spec():
    return np.arange(0, 1000, 0.01)  # density 100 per unit


def test_lambda_is_density_times_width(uniform_spec):
    t = Target("x", 500.0, 499.9, 500.1, "random")
    r = evaluate_target(uniform_spec, t)
    assert local_density(uniform_spec, 500.0) == pytest.approx(100, rel=1e-3)
    assert r["lam"] == pytest.approx(20, rel=1e-3)
    assert r["hits"] == pytest.approx(20, abs=1)


def test_miss_reports_nearest_distance(uniform_spec):
    t = Target("x", 500.005, 500.004, 500.006, "random")
    r = evaluate_target(uniform_spec, t)
    assert not r["hit"]
    assert r["best_rel"] == pytest.approx(0.005 / 500.005, rel=1e-6)


def test_correction_series_size_and_symmetry():
    c = correction_series(powers=range(7, 10))
    assert len(c) == 5**3
    np.testing.assert_allclose(c, -c[::-1])


def test_tuning_capacity_finds_shift():
    spec = np.array([10.0])
    corr = np.array([-0.5, 0.0, 0.25])
    assert tuning_capacity(spec, corr, Target("x", 10.25, 10.2499, 10.2501, "r"))
    assert not tuning_capacity(spec, corr, Target("x", 10.1, 10.0999, 10.1001, "r"))
