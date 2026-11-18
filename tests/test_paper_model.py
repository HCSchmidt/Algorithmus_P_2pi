"""
Checks against the paper itself (Paper.pdf, Fig. 3 code and worked examples).
"""

import math

import numpy as np
import pytest

from polynome2pi.energy_model import EnergyModel
from polynome2pi.evaluation.paper_model import (
    coefficient_grid,
    paper_energies,
    paper_energy_scalar,
)

# (tuple, value printed in / recomputed from the paper's Fig. 3 code, tolerance)
PAPER_EXAMPLES = {
    "proton": ((1, 1, 1, 1, 1, 1, 1), 1836.152673409, 1e-8),
    "neutron": ((1, 1, 1, 1, 1, 1, 0), 1838.68366245, 1e-8),
    "muon": ((0, 1, -1, 1, -1, 0, -1), 206.7682831, 1e-6),
    "pion_0_max": ((0, 1, 0, -2.5, -3, 0.5, 1), 264.1434, 1e-4),
    "pion_pm_max": ((0, 1, 1, 2.5, -3, 2.5, 0.5), 273.13259, 1e-5),
    "electron_1": ((0, 0, 0, 0.5, -2.5, 0.5, -0.5), 0.9957, 1e-4),
}


@pytest.mark.parametrize("name", PAPER_EXAMPLES)
def test_paper_formula_reproduces_paper_examples(name):
    coeffs, value, tol = PAPER_EXAMPLES[name]
    assert paper_energy_scalar(*coeffs) == pytest.approx(value, abs=tol)


def test_vectorised_paper_formula_matches_transcription():
    grid = coefficient_grid()
    idx = np.random.default_rng(0).choice(len(grid), 50_000, replace=False)
    sample = grid[idx]
    expected = np.array([paper_energy_scalar(*row) for row in sample])
    np.testing.assert_allclose(paper_energies(sample), expected, rtol=1e-14, atol=1e-11)


def test_paper_formula_base_parameter_replaces_pi_everywhere():
    row = np.array([[1, 1, 1, 1, 1, 1, 1]], dtype=float)
    b = 6.3
    assert paper_energies(row, base=b)[0] == pytest.approx(paper_energy_scalar(*row[0], pi=b / 2), rel=1e-14)


@pytest.mark.xfail(
    strict=True,
    reason=(
        "Known discrepancy: the repo's EnergyModel applies the E7 term to every "
        "zero/zero pair after pairing; the paper's Fig. 3 code applies it only to "
        "the first such pair per gluon, inside the loop. Remove this xfail once fixed."
    ),
)
@pytest.mark.parametrize("name", ["proton", "neutron", "muon"])
def test_repo_model_reproduces_paper_examples(name):
    coeffs, value, tol = PAPER_EXAMPLES[name]
    assert EnergyModel().energy(*coeffs) == pytest.approx(value, abs=tol)
