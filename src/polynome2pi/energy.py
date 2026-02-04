from __future__ import annotations

import cmath
from dataclasses import dataclass
from typing import Dict, Iterable, List, Tuple, Union

import numpy as np

from .cli import ScanSector
from .particles import Particle


PI = cmath.pi
TWO_PI = 2 * PI
TWO_PI_POW = {k: TWO_PI ** k for k in range(-30, 31)}

E_C_POS = (
    -PI
    + 2 * PI ** (-1)
    - PI ** (-3)
    + 2 * PI ** (-5)
    - PI ** (-7)
    + PI ** (-9)
    - PI ** (-12)
    - 2 * PI ** (-14)
)
E_C_NEG = 2 * PI - PI ** (-1) + E_C_POS
E_C_ZERO = PI ** (-12) + 2 * PI ** (-14)


@dataclass
class ParticleData:
    """
    Per-particle scan results keyed by m in [1..511].
    This mirrors the legacy arrays used by the original script so the CSV
    can stay identical.

    Conventions:
      - Emax/Emin store energies (float).
      - i_Emax/i_Emin store the scan index i_T where max/min occurred (int).
      - Dmax/Dmin store the raw integer coefficients (NOT divided by 2):
            [i4, i3, i2, i1, i0, i_minus_1, C]
        The report layer divides by 2 when writing CSV to match the legacy output.
    """

    Emax: np.ndarray      # shape (512,), float
    Emin: np.ndarray      # shape (512,), float
    i_Emax: np.ndarray    # shape (512,), int
    i_Emin: np.ndarray    # shape (512,), int
    Dmax: np.ndarray      # shape (512, 7), int (or float ok), raw loop ints
    Dmin: np.ndarray      # shape (512, 7), int (or float ok), raw loop ints

    @classmethod
    def empty(cls) -> "ParticleData":
        return cls(
            Emax=np.zeros(512, dtype=float),
            Emin=np.zeros(512, dtype=float),
            i_Emax=np.zeros(512, dtype=int),
            i_Emin=np.zeros(512, dtype=int),
            Dmax=np.zeros((512, 7), dtype=int),
            Dmin=np.zeros((512, 7), dtype=int),
        )

    def m_values(self) -> range:
        # legacy iterates m=1..511
        return range(1, 512)

class EnergieEngine:
    def __init__(self, sector: ScanSector):
        self.sector = sector
        self.E = [0.0] * 8
        self.g = [[0.0] * 10 for _ in range(10)]

        # map sector -> J4,J3,J2 (same as your presets)
        if sector is ScanSector.minimal:
            self.J4, self.J3, self.J2 = 0, 0, 1
        elif sector is ScanSector.light:
            self.J4, self.J3, self.J2 = 0, 1, 1
        elif sector is ScanSector.broad:
            self.J4, self.J3, self.J2 = 1, 1, 2
        elif sector is ScanSector.nucleon:
            self.J4, self.J3, self.J2 = 1, 1, 1
        elif sector is ScanSector.heavy:
            self.J4, self.J3, self.J2 = 2, 2, 2
        else:
            raise ValueError(f"Unhandled sector: {sector}")

    def energie(self, i4, i3, i2, i1, i0, i_1, C) -> float:
        g = self.g
        E = self.E

        g[2][4] = i4
        g[2][3] = i3
        g[2][2] = i2
        g[1][1] = i1
        g[1][0] = i0
        g[1][-1] = i_1

        if C > 0:
            E0 = C * E_C_POS
        elif C < 0:
            E0 = -C * E_C_NEG
        else:
            E0 = E_C_ZERO

        # gluonic
        for l in range(4, 1, -1):
            E0 += g[2][l] * TWO_PI_POW[l]

        # fermionic
        for n in range(1, -2, -1):
            E0 -= g[1][n] * TWO_PI_POW[n]

        # interaction terms
        for l in range(4, 1, -1):
            for n in range(1, -2, -1):
                if g[2][l] and g[1][n]:
                    if l + n < 4 and g[2][l] > 0:
                        E0 += g[2][l] * g[1][n] * 2 * TWO_PI_POW[-l - n - 1]
                    elif l + n < 4 and g[2][l] < 0:
                        E0 += g[2][l] * g[1][n] * 2 * TWO_PI_POW[-l - n]
                    elif l + n > 3:
                        E0 -= g[2][l] * g[1][n] * 2 * TWO_PI_POW[-l - n - 1]

                    E0 += abs(g[2][l] * g[1][n]) * 2 * TWO_PI_POW[-8]
                    g[2][l] = 0
                    g[1][n] = 0
                    break

        return float(E0)

    def _normalize_particles(
        self, particles: Union[Dict[str, Particle], Iterable[Particle]]
    ) -> List[Particle]:
        if isinstance(particles, dict):
            return list(particles.values())
        return list(particles)

    def run(self, particles: Union[Dict[str, Particle], Iterable[Particle]]):
        particles_list = self._normalize_particles(particles)

        xs_by_particle: Dict[str, List[int]] = {p.key: [] for p in particles_list}
        ys_by_particle: Dict[str, List[float]] = {p.key: [] for p in particles_list}
        grey_segments: List[Tuple[Tuple[int, float], Tuple[int, float]]] = []

        # per-particle storage for report
        results: Dict[str, ParticleData] = {}
        for p in particles_list:
            results[p.key] = ParticleData(
                Emax=np.zeros(512, dtype=float),
                Emin=np.zeros(512, dtype=float),
                i_Emax=np.zeros(512, dtype=int),
                i_Emin=np.zeros(512, dtype=int),
                Dmax=np.zeros((512, 7), dtype=float),
                Dmin=np.zeros((512, 7), dtype=float),
            )

        # counts per m (legacy logic)
        Cnt = [0] * 520
        mmax = 0
        ct = 0

        labels: List[Tuple[ScanSector, Particle, int, float, str]] = []
        labeled_blocks = set()

        i_T = 0
        i_T1 = 0

        for i4 in range(-2 * self.J4, 2 * self.J4 + 1):
            for i3 in range(-2 * self.J3, 2 * self.J3 + 1):
                for i2 in range(-2 * self.J2, 2 * self.J2 + 1):
                    for i1 in range(-6, 7):
                        for i0 in range(-6, 7):
                            for i_1 in range(-6, 7):
                                for C in range(-2, 3):
                                    E0 = self.energie(
                                        i4 / 2, i3 / 2, i2 / 2,
                                        i1 / 2, i0 / 2, i_1 / 2,
                                        C / 2
                                    )

                                    if E0 <= 0:
                                        continue

                                    if self.sector is ScanSector.heavy and E0 < 1500:
                                        continue
                                    if self.sector is ScanSector.nucleon and not (1836 <= E0 <= 1839):
                                        continue

                                    i_T += 1

                                    m = int(256 + 32 * i4 + 4 * i3 + i2)
                                    if m > mmax:
                                        mmax = m
                                        ct = 0
                                    ct += 1
                                    Cnt[mmax] = ct

                                    matched_any = False

                                    for p in particles_list:
                                        d = E0 - p.theory_E
                                        if d <= p.sd_plus and d >= p.sd_minus:
                                            matched_any = True
                                            i_T1 += 1

                                            pdata = results[p.key]

                                            if pdata.Emax[m] <= E0:
                                                pdata.Emax[m] = E0
                                                pdata.i_Emax[m] = i_T
                                                pdata.Dmax[m] = [i4, i3, i2, i1, i0, i_1, C]

                                            if pdata.Emin[m] == 0 or pdata.Emin[m] >= E0:
                                                pdata.Emin[m] = E0
                                                pdata.i_Emin[m] = i_T
                                                pdata.Dmin[m] = [i4, i3, i2, i1, i0, i_1, C]

                                            xs_by_particle[p.key].append(i_T)
                                            ys_by_particle[p.key].append(E0)

                                            # add one label per (particle, m) block (like old "flag==1")
                                            key = (p.key, m)
                                            if key not in labeled_blocks:
                                                labeled_blocks.add(key)
                                                e_max = pdata.Emax[m]
                                                e_min = pdata.Emin[m] if pdata.Emin[m] != 0 else e_max
                                                mean = (e_max + e_min) / 2
                                                labels.append((self.sector, p, pdata.i_Emin[m], mean))

                                    if not matched_any:
                                        grey_segments.append(((i_T, E0), (i_T + 1, E0)))

        return (
            xs_by_particle,
            ys_by_particle,
            grey_segments,
            results,
            Cnt,
            labels,
            i_T,
            i_T1,
        )