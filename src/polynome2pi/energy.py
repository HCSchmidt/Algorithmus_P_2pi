from __future__ import annotations

import cmath
from dataclasses import dataclass
from typing import Dict, List, Tuple

import numpy as np

from .cli import ScanSector
from .particles import Particle


# ---------------------------------------------------------------------
# Internal constants (precomputed once)
# ---------------------------------------------------------------------

PI = cmath.pi
TWO_PI = 2 * PI

# powers used in Energie()
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


# ---------------------------------------------------------------------
# Result container
# ---------------------------------------------------------------------

@dataclass
class ScanResult:
    particle: Particle
    m_index: int
    E_min: float
    E_max: float
    i_min: int
    i_max: int
    D: Tuple[float, float, float, float, float, float, float]


# ---------------------------------------------------------------------
# Energy engine
# ---------------------------------------------------------------------

class EnergieEngine:
    """
    Stateless energy scan engine.
    """

    def __init__(self, sector: ScanSector):
        self.sector = sector

        # working buffers (local to instance, no globals)
        self.E = [0.0] * 8
        self.g = [[0.0] * 10 for _ in range(10)]

    # -----------------------------------------------------------------
    # Core energy function (unchanged physics)
    # -----------------------------------------------------------------
    def energie(
        self,
        i4, i3, i2, i1, i0, i_1, C
    ) -> float:
        g = self.g
        E = self.E

        g[2][4] = i4
        g[2][3] = i3
        g[2][2] = i2
        g[1][1] = i1
        g[1][0] = i0
        g[1][-1] = i_1

        # C-dependent base term
        if C > 0:
            E[0] = C * E_C_POS
        elif C < 0:
            E[0] = -C * E_C_NEG
        else:
            E[0] = E_C_ZERO

        # gluonic
        for l in range(4, 1, -1):
            E[0] += g[2][l] * TWO_PI_POW[l]

        # fermionic
        for n in range(1, -2, -1):
            E[0] -= g[1][n] * TWO_PI_POW[n]

        # interaction terms
        for l in range(4, 1, -1):
            for n in range(1, -2, -1):
                if g[2][l] and g[1][n]:
                    if l + n < 4 and g[2][l] > 0:
                        E[0] += g[2][l] * g[1][n] * 2 * TWO_PI_POW[-l - n - 1]
                    elif l + n < 4 and g[2][l] < 0:
                        E[0] += g[2][l] * g[1][n] * 2 * TWO_PI_POW[-l - n]
                    elif l + n > 3:
                        E[0] -= g[2][l] * g[1][n] * 2 * TWO_PI_POW[-l - n - 1]

                    E[0] += abs(g[2][l] * g[1][n]) * 2 * TWO_PI_POW[-8]
                    g[2][l] = 0
                    g[1][n] = 0
                    break

        return float(E[0])

    # -----------------------------------------------------------------
    # Main scan
    # -----------------------------------------------------------------
    def run(
        self,
        particles: List[Particle],
    ):
        """
        Run full scan.
        Returns everything needed for plotting and reporting.
        """

        xs_by_particle: Dict[str, List[int]] = {p.key: [] for p in particles}
        ys_by_particle: Dict[str, List[float]] = {p.key: [] for p in particles}
        grey_segments: List[Tuple[Tuple[int, float], Tuple[int, float]]] = []

        results: Dict[str, List[ScanResult]] = {p.key: [] for p in particles}
        counts: Dict[int, int] = {}

        i_T = 0
        i_T1 = 0

        for i4 in range(-2, 3):
            for i3 in range(-2, 3):
                for i2 in range(-2, 3):
                    for i1 in range(-6, 7):
                        for i0 in range(-6, 7):
                            for i_1 in range(-6, 7):
                                for C in range(-2, 3):
                                    i_T += 1

                                    E0 = self.energie(
                                        i4 / 2, i3 / 2, i2 / 2,
                                        i1 / 2, i0 / 2, i_1 / 2,
                                        C / 2,
                                    )

                                    if E0 <= 0:
                                        continue

                                    # sector cuts
                                    if self.sector is ScanSector.heavy and E0 < 1500:
                                        continue
                                    if self.sector is ScanSector.nucleon and not (1836 <= E0 <= 1839):
                                        continue

                                    i_T1 += 1

                                    matched = False
                                    for p in particles:
                                        if p.matches(E0):
                                            xs_by_particle[p.key].append(i_T)
                                            ys_by_particle[p.key].append(E0)
                                            matched = True

                                    if not matched:
                                        grey_segments.append(((i_T, E0), (i_T + 1, E0)))

        return (
            xs_by_particle,
            ys_by_particle,
            grey_segments,
            results,
            counts,
            i_T,
            i_T1,
        )