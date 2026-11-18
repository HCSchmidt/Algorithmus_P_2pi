from __future__ import annotations

import cmath
from dataclasses import dataclass


@dataclass(frozen=True)
class EnergyModel:
    """
    Implements the legacy energy computation E(i4,i3,i2,i1,i0,i-1,C).

    Sensitivity analysis: vary ONLY the base (2π) by a scale factor,
    while keeping π itself constant (so π^-k terms remain unchanged).
    """

    base_scale: float = 1.0

    @property
    def pi(self) -> complex:
        return cmath.pi

    @property
    def two_pi(self) -> complex:
        # only (2π) is scaled
        return (2 * self.pi) * float(self.base_scale)

    def energy(self, i4, i3, i2, i1, i0, i_minus1, C) -> float:
        pi = self.pi
        two_pi = self.two_pi

        # --- C-dependent constants (KEEP pi unscaled) ---
        E_C_pos = (
            -pi
            + 2 * pi ** (-1)
            - pi ** (-3)
            + 2 * pi ** (-5)
            - pi ** (-7)
            + pi ** (-9)
            - pi ** (-12)
            - 2 * pi ** (-14)
        )
        E_C_neg = 2 * pi - pi ** (-1) + E_C_pos
        E0 = pi ** (-12) + 2 * pi ** (-14)

        if C > 0:
            E0 = C * E_C_pos
        elif C < 0:
            E0 = -C * E_C_neg

        # --- Gluon-like part (uses scaled 2π) ---
        E2 = i4 * (two_pi**4) + i3 * (two_pi**3) + i2 * (two_pi**2)

        # --- Fermion-like part (uses scaled 2π) ---
        E1 = -(i1 * (two_pi**1) + i0 * (two_pi**0) + i_minus1 * (two_pi ** (-1)))

        # --- Interaction terms (uses scaled 2π) ---
        E3 = 0.0
        E4 = 0.0
        E5 = 0.0
        E6 = 0.0
        E7 = 0.0

        gluons = [(4, i4), (3, i3), (2, i2)]
        fermions = [(1, i1), (0, i0), (-1, i_minus1)]

        g = [(l, val) for l, val in gluons]
        f = [(n, val) for n, val in fermions]

        for gi, (l, gv) in enumerate(g):
            if gv == 0:
                continue
            for fi, (n, fv) in enumerate(f):
                if fv == 0:
                    continue

                if (l + n) < 4 and gv > 0:
                    E3 += gv * fv * 2 * (two_pi ** (-l - n - 1))
                if (l + n) < 4 and gv < 0:
                    E4 += gv * fv * 2 * (two_pi ** (-l - n))
                if (l + n) > 3:
                    E5 -= gv * fv * 2 * (two_pi ** (-l - n - 1))

                E6 += abs(gv * fv) * 2 * (two_pi ** (-8))

                g[gi] = (l, 0.0)
                f[fi] = (n, 0.0)
                break

        for l, gv in g:
            for n, fv in f:
                if gv == 0 and fv == 0:
                    E7 -= two_pi ** (-l - n - 1)
                    E7 -= two_pi ** (-l - n)

        E = E0 + E1 + E2 + E3 + E4 + E5 + E6 + E7
        return float(E.real)
