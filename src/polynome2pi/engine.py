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
    possible_ET: int
    real_ET: int


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
        model = self.model

        # bins: particle_key -> m -> BinResult under construction
        bins: Dict[str, Dict[int, BinResult]] = {k: {} for k in particles.keys()}

        # counts per m (legacy used Cnt[mmax], here: counts per m directly)
        m_counts = np.zeros(520, dtype=int)

        # plotting buffers
        matched_xs: Dict[str, List[int]] = {k: [] for k in particles.keys()}
        matched_ys: Dict[str, List[float]] = {k: [] for k in particles.keys()}
        unmatched_segments: List[Tuple[Tuple[int, float], Tuple[int, float]]] = []

        possible_ET = 0
        real_ET = 0
        i_T = 0  # scan index over accepted E0 >=0 combinations

        # hot locals
        energy = model.energy
        m_index = self._m_index

        # scan loop (mirrors legacy ranges)
        for i4 in range(-2 * preset.J4, 2 * preset.J4 + 1):
            for i3 in range(-2 * preset.J3, 2 * preset.J3 + 1):
                for i2 in range(-2 * preset.J2, 2 * preset.J2 + 1):
                    for i1 in range(-6, 7):
                        for i0 in range(-6, 7):
                            for i_minus1 in range(-6, 7):
                                for C in range(-2, 3):
                                    # legacy uses halves
                                    E0 = energy(
                                        i4 / 2,
                                        i3 / 2,
                                        i2 / 2,
                                        i1 / 2,
                                        i0 / 2,
                                        i_minus1 / 2,
                                        C / 2,
                                    )
                                    possible_ET += 1
                                    if E0 < 0:
                                        continue

                                    # sector-specific cuts (legacy-like)
                                    if preset.sector.value == "heavy" and E0 < 1500:
                                        continue
                                    if preset.sector.value == "nucleon" and (E0 < 1836 or E0 > 1839):
                                        continue

                                    i_T += 1
                                    m = m_index(i4, i3, i2)
                                    if 0 <= m < len(m_counts):
                                        m_counts[m] += 1

                                    # match against particles
                                    matched_any = False
                                    for pkey, p in particles.items():
                                        # tolerance window around theory_E
                                        if (E0 - p.theory_E) <= p.sd_plus and (E0 - p.theory_E) >= p.sd_minus:
                                            matched_any = True
                                            real_ET += 1

                                            # update plotting buffers
                                            matched_xs[pkey].append(i_T)
                                            matched_ys[pkey].append(E0)

                                            # update per-particle, per-m-bin results
                                            coeffs: Coeff7 = (i4, i3, i2, i1, i0, i_minus1, C)
                                            current = bins[pkey].get(m)
                                            if current is None:
                                                bins[pkey][m] = BinResult(
                                                    m=m,
                                                    i_T_min=i_T,
                                                    i_T_max=i_T,
                                                    E_min=E0,
                                                    E_max=E0,
                                                    coeff_min=coeffs,
                                                    coeff_max=coeffs,
                                                    counts=int(m_counts[m]),
                                                )
                                            else:
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
                                                current.counts = int(m_counts[m])

                                    if not matched_any:
                                        unmatched_segments.append(((i_T, E0), (i_T + 1, E0)))

        bins_by_particle: Dict[str, List[BinResult]] = {
            pkey: sorted(mmap.values(), key=lambda r: r.m) for pkey, mmap in bins.items()
        }

        return ScanOutputs(
            matched_points={k: (matched_xs[k], matched_ys[k]) for k in particles.keys()},
            unmatched_segments=unmatched_segments,
            bins_by_particle=bins_by_particle,
            possible_ET=i_T,  # comparable to legacy “possible ET” after cuts
            real_ET=real_ET,
        )