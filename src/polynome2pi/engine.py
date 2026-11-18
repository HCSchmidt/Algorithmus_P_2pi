from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Tuple

import numpy as np

from .particles import Particle
from .presets import ScanPreset
from .energy_model import EnergyModel


Coeff7 = Tuple[int, int, int, int, int, int, int]  # i4,i3,i2,i1,i0,i-1,C


@dataclass
class BinResult:
    """
    One 'm-bin' match segment for a particle:
      - where we saw it first/last (i_T indices)
      - min/max energy within that bin
      - coefficients for min/max
      - number of times this bin occurred (counts)
    """

    m: int
    i_T_min: int
    i_T_max: int
    E_min: float
    E_max: float
    coeff_min: Coeff7
    coeff_max: Coeff7
    counts: int

    @property
    def delta_i(self) -> int:
        return self.i_T_max - self.i_T_min + 1

    @property
    def mean_E(self) -> float:
        return (self.E_min + self.E_max) / 2.0

    @property
    def delta_i_over_counts(self) -> float:
        return (self.delta_i / self.counts) if self.counts else 0.0


@dataclass
class ScanOutputs:
    # For plotting (batched)
    matched_points: Dict[str, Tuple[List[int], List[float]]]  # particle_key -> (xs, ys)
    unmatched_segments: List[Tuple[Tuple[int, float], Tuple[int, float]]]  # grey segments

    # Results
    bins_by_particle: Dict[str, List[BinResult]]

    # Global scan stats
    accepted_scan_points: int
    total_particle_hits: int

    
    grid_shape: Tuple[int, int, int]                 # (n_i4, n_i3, n_i2)
    grid_i4_values: List[int]                        # actual i4 integer values scanned
    grid_i3_values: List[int]
    grid_i2_values: List[int]
    any_match_grid: np.ndarray                       # shape (n_i4, n_i3, n_i2) counts of "matched_any"
    particle_match_grids: Dict[str, np.ndarray]      # per particle matched count grid


@dataclass
class BinStore:
    """Small wrapper around per-particle m-bin results with clear intent."""

    _by_m: Dict[int, BinResult]

    def get(self, m: int) -> BinResult | None:
        return self._by_m.get(m)

    def upsert(self, m: int, i_T: int, E0: float, coeffs: Coeff7, m_count: int) -> None:
        current = self._by_m.get(m)
        if current is None:
            self._by_m[m] = BinResult(
                m=m,
                i_T_min=i_T,
                i_T_max=i_T,
                E_min=E0,
                E_max=E0,
                coeff_min=coeffs,
                coeff_max=coeffs,
                counts=int(m_count),
            )
            return

        # keep first/last occurrence in i_T
        current.i_T_min = min(current.i_T_min, i_T)
        current.i_T_max = max(current.i_T_max, i_T)

        # update min/max energy and coeffs
        if E0 <= current.E_min or current.E_min == 0:
            current.E_min = E0
            current.coeff_min = coeffs
        if E0 >= current.E_max:
            current.E_max = E0
            current.coeff_max = coeffs

        # update counts
        current.counts = int(m_count)

    def values_sorted(self) -> List[BinResult]:
        return sorted(self._by_m.values(), key=lambda r: r.m)


@dataclass
class ParticleAccumulator:
    """All scan outputs we collect for one particle."""

    particle: Particle
    bins: BinStore
    xs: List[int]
    ys: List[float]

    def record_match(self, i_T: int, m: int, E0: float, coeffs: Coeff7, m_count: int) -> None:
        # Plot buffers
        self.xs.append(i_T)
        self.ys.append(E0)

        # Result bins
        self.bins.upsert(m=m, i_T=i_T, E0=E0, coeffs=coeffs, m_count=m_count)


class ScanEngine:
    def __init__(self, preset: ScanPreset, model: EnergyModel):
        self.preset = preset
        self.model = model

    @staticmethod
    def _m_index(i4: int, i3: int, i2: int) -> int:
        # legacy mapping
        return int(256 + 32 * i4 + 4 * i3 + i2)

    def run(self, particles: Dict[str, Particle]) -> ScanOutputs:
        preset = self.preset

        i4_values = list(range(-2 * preset.J4, 2 * preset.J4 + 1))
        i3_values = list(range(-2 * preset.J3, 2 * preset.J3 + 1))
        i2_values = list(range(-2 * preset.J2, 2 * preset.J2 + 1))
        grid_shape = (len(i4_values), len(i3_values), len(i2_values))

        any_match_grid = np.zeros(grid_shape, dtype=int)
        particle_match_grids = {k: np.zeros(grid_shape, dtype=int) for k in particles}
        accumulators: Dict[str, ParticleAccumulator] = {
            key: ParticleAccumulator(particle=p, bins=BinStore(_by_m={}), xs=[], ys=[])
            for key, p in particles.items()
        }
        m_counts = np.zeros(520, dtype=int)
        unmatched_segments: List[Tuple[Tuple[int, float], Tuple[int, float]]] = []
        total_particle_hits = 0
        scan_index = 0

        # One chunk per i4 value keeps memory bounded and preserves scan order.
        for gi4, i4 in enumerate(i4_values):
            chunk = self._scan_chunk(i4, i3_values, i2_values)
            if chunk is None:
                continue
            E, coeffs, grid_idx, m = chunk
            n = len(E)
            scan_ids = np.arange(scan_index + 1, scan_index + n + 1)
            scan_index += n

            m_count_at_point = self._running_counts(m, m_counts)

            matched_any = np.zeros(n, dtype=bool)
            for pkey, acc in accumulators.items():
                p = acc.particle
                diff = E - p.theory_E
                hit = (diff <= p.sd_plus) & (diff >= p.sd_minus)
                idx = np.flatnonzero(hit)
                if idx.size == 0:
                    continue
                matched_any |= hit
                total_particle_hits += int(idx.size)
                np.add.at(particle_match_grids[pkey], (gi4, grid_idx[0][idx], grid_idx[1][idx]), 1)
                for i_T, mm, e, c, cnt in zip(
                    scan_ids[idx].tolist(),
                    m[idx].tolist(),
                    E[idx].tolist(),
                    coeffs[idx].tolist(),
                    m_count_at_point[idx].tolist(),
                ):
                    acc.record_match(i_T=i_T, m=mm, E0=e, coeffs=tuple(c), m_count=cnt)

            hit_idx = np.flatnonzero(matched_any)
            np.add.at(any_match_grid, (gi4, grid_idx[0][hit_idx], grid_idx[1][hit_idx]), 1)

            miss_idx = np.flatnonzero(~matched_any)
            unmatched_segments.extend(
                ((s, e), (s + 1, e))
                for s, e in zip(scan_ids[miss_idx].tolist(), E[miss_idx].tolist())
            )

        return ScanOutputs(
            matched_points={key: (acc.xs, acc.ys) for key, acc in accumulators.items()},
            unmatched_segments=unmatched_segments,
            bins_by_particle={key: acc.bins.values_sorted() for key, acc in accumulators.items()},
            accepted_scan_points=scan_index,
            total_particle_hits=total_particle_hits,
            grid_shape=grid_shape,
            grid_i4_values=i4_values,
            grid_i3_values=i3_values,
            grid_i2_values=i2_values,
            any_match_grid=any_match_grid,
            particle_match_grids=particle_match_grids,
        )

    def _scan_chunk(self, i4: int, i3_values: List[int], i2_values: List[int]):
        """Energies of all accepted points with this i4, in legacy loop order."""
        f = np.arange(-6, 7)
        c = np.arange(-2, 3)
        g3, g2, g1, g0, gm1, gc = (
            a.ravel()
            for a in np.meshgrid(
                np.arange(len(i3_values)), np.arange(len(i2_values)), f, f, f, c, indexing="ij"
            )
        )
        i3 = np.asarray(i3_values)[g3]
        i2 = np.asarray(i2_values)[g2]
        i4_arr = np.full_like(i3, i4)

        E = self.model.energies(i4_arr / 2, i3 / 2, i2 / 2, g1 / 2, g0 / 2, gm1 / 2, gc / 2)

        accepted = ~(E < 0)
        # sector-specific cuts (legacy-like); note these compare against preset
        # names, not ScanSector values, so they never trigger with the current enum.
        sector = self.preset.sector.value
        if sector == "heavy":
            accepted &= ~(E < 1500)
        if sector == "nucleon":
            accepted &= ~((E < 1836) | (E > 1839))

        idx = np.flatnonzero(accepted)
        if idx.size == 0:
            return None
        coeffs = np.stack([i4_arr, i3, i2, g1, g0, gm1, gc], axis=1)[idx]
        m = 256 + 32 * i4 + 4 * i3[idx] + i2[idx]
        return E[idx], coeffs, (g3[idx], g2[idx]), m

    @staticmethod
    def _running_counts(m: np.ndarray, m_counts: np.ndarray) -> np.ndarray:
        """Per point: how often its m has been seen so far (inclusive). Updates m_counts."""
        if m.min() < 0 or m.max() >= len(m_counts):
            raise ValueError("m index out of legacy range")
        order = np.argsort(m, kind="stable")
        sorted_m = m[order]
        group_start = np.r_[0, np.flatnonzero(np.diff(sorted_m)) + 1]
        rank = np.arange(len(m)) - np.repeat(group_start, np.diff(np.r_[group_start, len(m)]))
        counts = np.empty_like(m)
        counts[order] = m_counts[sorted_m] + rank + 1
        np.add.at(m_counts, m, 1)
        return counts
