from __future__ import annotations

import cmath
from dataclasses import dataclass


@dataclass(frozen=True)
class EnergyModel:
    """
    Implements the legacy energy computation E(i4,i3,i2,i1,i0,i-1,C)
    with precomputed constants.
    """

    pi: complex = cmath.pi

    def __post_init__(self):
        # dataclass frozen: use object.__setattr__ if needed (not needed here)
        pass

    @property
    def two_pi(self) -> complex:
        return 2 * self.pi

    def energy(self, i4, i3, i2, i1, i0, i_minus1, C) -> float:
        pi = self.pi
        two_pi = 2 * pi

        # Precomputed-like constants (kept inline for clarity; fast enough)
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

        # Gluon-like part (l=4..2)
        E2 = i4 * (two_pi ** 4) + i3 * (two_pi ** 3) + i2 * (two_pi ** 2)

        # Fermion-like part (n=1..-1)
        E1 = -(i1 * (two_pi ** 1) + i0 * (two_pi ** 0) + i_minus1 * (two_pi ** (-1)))

        # Interaction terms replicate the legacy logic: pair cancellations
        # We implement the same structure but without mutating shared arrays.
        E3 = 0.0  # neutral matter
        E4 = 0.0  # neutral antimatter
        E5 = 0.0  # gravitation-like
        E6 = 0.0  # internal time
        E7 = 0.0  # additional neutral terms

        # Legacy code “pairs” non-zero gluon coefficients with fermion coefficients.
        gluons = [(4, i4), (3, i3), (2, i2)]
        fermions = [(1, i1), (0, i0), (-1, i_minus1)]

        # Convert to mutable copies to emulate "consume once then break"
        g = [(l, val) for l, val in gluons]
        f = [(n, val) for n, val in fermions]

        for gi, (l, gv) in enumerate(g):
            if gv == 0:
                continue
            for fi, (n, fv) in enumerate(f):
                if fv == 0:
                    continue

                # replicate boolean multipliers (cast to int via Python truthiness)
                if (l + n) < 4 and gv > 0:
                    E3 += gv * fv * 2 * (two_pi ** (-l - n - 1))
                if (l + n) < 4 and gv < 0:
                    E4 += gv * fv * 2 * (two_pi ** (-l - n))
                if (l + n) > 3:
                    E5 -= gv * fv * 2 * (two_pi ** (-l - n - 1))

                E6 += abs(gv * fv) * 2 * (two_pi ** (-8))

                # consume both once (legacy sets to 0 and breaks)
                g[gi] = (l, 0.0)
                f[fi] = (n, 0.0)
                break

        # If both are zero in some combinations, legacy subtracts extra neutral terms;
        # we approximate that behaviour for remaining zero pairs.
        for (l, gv) in g:
            for (n, fv) in f:
                if gv == 0 and fv == 0:
                    E7 -= (two_pi ** (-l - n - 1))
                    E7 -= (two_pi ** (-l - n))

        E = E0 + E1 + E2 + E3 + E4 + E5 + E6 + E7

        # Return real part as float (legacy uses cmath but energies are treated real)
        return float(E.real)