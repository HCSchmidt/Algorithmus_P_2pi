"""
Null-hypothesis evaluation of the P(2π) mass formula.

The spectrum of a formula is the sorted set of all energies over the full
coefficient grid (union of all sectors). For a target window [lo, hi]:

  hit      : some energy lies in the window
  lam      : expected number of energies in the window if the spectrum were
             locally uniform = local density × window width (chance level)
  best_rel : relative distance of the nearest energy to the target
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict, Iterable, List

import numpy as np

from polynome2pi.evaluation.paper_model import coefficient_grid, paper_energies
from polynome2pi.particles import get_particles

M_E_MEV = 0.51099895000


@dataclass(frozen=True)
class Target:
    key: str
    mass: float  # m_e
    lo: float
    hi: float
    group: str  # "table" | "holdout"


def table_targets() -> List[Target]:
    return [
        Target(k, p.theory_E, p.theory_E + p.sd_minus, p.theory_E + p.sd_plus, "table")
        for k, p in get_particles().items()
    ]


# PDG masses (MeV, 1σ) of well-measured hadrons that are NOT in the repo's table,
# plus positronium (2 m_e minus 6.8029 eV binding). Used as an out-of-sample test.
_HOLDOUT_MEV = {
    "lambda": (1115.683, 0.006),
    "sigma_plus": (1189.37, 0.07),
    "sigma_0": (1192.642, 0.024),
    "sigma_minus": (1197.449, 0.030),
    "xi_0": (1314.86, 0.20),
    "xi_minus": (1321.71, 0.07),
    "sigma_star_plus": (1382.83, 0.34),
    "xi_star_0": (1531.80, 0.32),
    "omega_minus": (1672.45, 0.29),
}


def holdout_targets() -> List[Target]:
    out = [
        Target(k, m / M_E_MEV, (m - s) / M_E_MEV, (m + s) / M_E_MEV, "holdout")
        for k, (m, s) in _HOLDOUT_MEV.items()
    ]
    # Bound states known far more precisely than the scan's level spacing.
    # Binding energies in eV -> m_e; constituent uncertainties propagated.
    ev = 1e-6 / M_E_MEV
    mp, sp = 1836.152673426, 3.2e-8
    mmu, smu = 206.7682827, 4.6e-6
    md, sd = 3670.4829677, 1.1e-6
    sharp = {
        "positronium": (2 - 6.8029 * ev, 1e-9),
        "hydrogen_1H": (mp + 1 - 13.5984 * ev, sp),
        "muonium": (mmu + 1 - 13.5390 * ev, smu),
        "deuterium_atom": (md + 1 - 13.6021 * ev, sd),
    }
    out += [Target(k, m, m - s, m + s, "holdout") for k, (m, s) in sharp.items()]
    return out


def pdg_electron() -> Target:
    """Electron with its real precision: exactly 1 m_e by definition."""
    return Target("electron_exact", 1.0, 1 - 1e-9, 1 + 1e-9, "pdg")


def spectrum(base: float = 2 * math.pi, grid: np.ndarray | None = None) -> np.ndarray:
    g = coefficient_grid() if grid is None else grid
    E = paper_energies(g, base=base)
    return np.sort(E[E >= 0])


def local_density(spec: np.ndarray, x: float) -> float:
    half = max(0.02 * x, 0.5)
    n = np.searchsorted(spec, x + half) - np.searchsorted(spec, x - half)
    return n / (2 * half)


def evaluate_target(spec: np.ndarray, t: Target) -> Dict:
    n_in = int(np.searchsorted(spec, t.hi, "right") - np.searchsorted(spec, t.lo, "left"))
    j = np.searchsorted(spec, t.mass)
    near = spec[max(j - 1, 0) : j + 1]
    best = float(near[np.argmin(np.abs(near - t.mass))])
    in_range = spec[0] <= t.mass <= spec[-1]
    lam = local_density(spec, t.mass) * (t.hi - t.lo)
    return {
        "key": t.key,
        "group": t.group,
        "mass": t.mass,
        "width": t.hi - t.lo,
        "in_range": bool(in_range),
        "hits": n_in,
        "hit": n_in > 0,
        "lam": lam,
        "p_chance": 1 - math.exp(-lam),
        "best": best,
        "best_rel": abs(best - t.mass) / t.mass,
    }


def evaluate_targets(spec: np.ndarray, targets: Iterable[Target]) -> List[Dict]:
    return [evaluate_target(spec, t) for t in targets]


def base_null(
    targets: List[Target], bases: np.ndarray, grid: np.ndarray | None = None
) -> List[Dict]:
    """For every base: which targets are hit, and their nearest relative distance."""
    g = coefficient_grid() if grid is None else grid
    rows = []
    for b in bases:
        spec = spectrum(float(b), g)
        for r in evaluate_targets(spec, targets):
            rows.append({"base": float(b), **r})
    return rows


def correction_series(pi: float = math.pi, powers=range(7, 15), coeffs=(-2, -1, 0, 1, 2)) -> np.ndarray:
    """Sorted values of Σ a_k π^-k — the same kind of term the paper's E_C constants are built from."""
    vals = np.zeros(1)
    for k in powers:
        vals = (vals[:, None] + np.array(coeffs)[None, :] * pi ** (-k)).ravel()
    return np.sort(vals)


def tuning_capacity(spec: np.ndarray, corr: np.ndarray, t: Target, neighbours: int = 50) -> bool:
    """
    Can the target be hit by a grid energy plus ONE extra correction constant
    from `corr`? Mirrors how E_C_pos / E0 shift a grid value by Σ a_k π^-k.
    """
    j = np.searchsorted(spec, t.mass)
    for E in spec[max(j - neighbours, 0) : j + neighbours]:
        lo, hi = t.lo - E, t.hi - E
        if np.searchsorted(corr, hi, "right") > np.searchsorted(corr, lo, "left"):
            return True
    return False


def random_targets(
    spec: np.ndarray, rel_widths: np.ndarray, n: int, lo: float, hi: float, seed: int = 0
) -> List[Target]:
    """Log-uniform random masses in [lo, hi], each with a relative width drawn from `rel_widths`."""
    rng = np.random.default_rng(seed)
    masses = np.exp(rng.uniform(math.log(lo), math.log(hi), n))
    w = rng.choice(rel_widths, n) * masses
    return [Target(f"rand{i}", m, m - x / 2, m + x / 2, "random") for i, (m, x) in enumerate(zip(masses, w))]
