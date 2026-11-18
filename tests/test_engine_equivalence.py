"""
The vectorised EnergyModel/ScanEngine must reproduce the original
implementation (frozen in tests/reference_*.py) exactly — not approximately.
"""

import filecmp
import itertools
from dataclasses import asdict
from pathlib import Path

import numpy as np
import pytest
import reference_energy_model
import reference_engine

from polynome2pi.cli import ScanSector
from polynome2pi.energy_model import EnergyModel
from polynome2pi.engine import ScanEngine
from polynome2pi.particles import get_particles
from polynome2pi.presets import preset_for_sector
from polynome2pi.scan.run import run_scan

REPO = Path(__file__).resolve().parents[1]


def _random_coeffs(n, seed=0):
    rng = np.random.default_rng(seed)
    g = rng.integers(-4, 5, size=(n, 3))
    f = rng.integers(-6, 7, size=(n, 3))
    c = rng.integers(-2, 3, size=(n, 1))
    return np.hstack([g, f, c]) / 2


@pytest.mark.parametrize("base_scale", [1.0, 0.999, 1.0013])
def test_vectorised_energy_is_bit_identical(base_scale):
    coeffs = _random_coeffs(50_000)
    # include every zero-pattern explicitly: pairing logic branches on it
    patterns = np.array(list(itertools.product([0, 1, -1], repeat=6)))
    patterns = np.hstack([patterns / 2, np.zeros((len(patterns), 1))])
    coeffs = np.vstack([coeffs, patterns])

    ref = reference_energy_model.EnergyModel(base_scale=base_scale)
    new = EnergyModel(base_scale=base_scale)
    expected = np.array([ref.energy(*row) for row in coeffs])
    actual = new.energies(*coeffs.T)
    assert np.array_equal(actual, expected)


def _assert_outputs_equal(a, b):
    assert a.accepted_scan_points == b.accepted_scan_points
    assert a.total_particle_hits == b.total_particle_hits
    assert a.grid_shape == b.grid_shape
    assert a.unmatched_segments == b.unmatched_segments
    assert a.matched_points == b.matched_points
    # BinResult is a different class in the frozen reference, so compare field values
    as_dicts = lambda bins: {k: [asdict(r) for r in v] for k, v in bins.items()}
    assert as_dicts(a.bins_by_particle) == as_dicts(b.bins_by_particle)
    assert np.array_equal(a.any_match_grid, b.any_match_grid)
    for k in a.particle_match_grids:
        assert np.array_equal(a.particle_match_grids[k], b.particle_match_grids[k])


@pytest.mark.parametrize("sector", [ScanSector.MINIMAL, ScanSector.LIGHT])
def test_engine_matches_reference(sector):
    preset = preset_for_sector(sector)
    particles = get_particles()
    expected = reference_engine.ScanEngine(preset, reference_engine.EnergyModel()).run(particles)
    actual = ScanEngine(preset, EnergyModel()).run(particles)
    _assert_outputs_equal(actual, expected)


@pytest.mark.parametrize("sector", [ScanSector.MINIMAL, ScanSector.LIGHT, ScanSector.BROAD])
def test_committed_results_reproduce(sector, tmp_path):
    run_scan(sector, None, tmp_path)
    committed = REPO / "results" / sector.value
    for ext in (".csv", ".txt"):
        name = f"scan_{sector.value}{ext}"
        assert filecmp.cmp(tmp_path / name, committed / name, shallow=False), name
