"""
The energy formula exactly as printed in the paper (Fig. 3, `Energie(...)`),
vectorised and with a configurable base.

Differences to `polynome2pi.energy_model.EnergyModel` (the repo port):
the paper applies the E7 term inside the pairing loop, only for the first
(l, n) pair where both coefficients are zero, followed by `break`. The repo
port subtracts E7 for every zero/zero pair after pairing, which shifts most
energies (e.g. the proton tuple gives 1835.894 instead of 1836.1527).
"""

from __future__ import annotations

import math

import numpy as np

GLUON_POWERS = (4, 3, 2)
FERMION_POWERS = (1, 0, -1)


def paper_energy_scalar(i4, i3, i2, i1, i0, i_minus1, C, pi: float = math.pi) -> float:
    """Literal transcription of the paper's Fig. 3 code (reference for tests)."""
    g2 = {4: i4, 3: i3, 2: i2}
    g1 = {1: i1, 0: i0, -1: i_minus1}
    E = [0.0] * 8
    E_C_pos = (
        -pi + 2 * pi**-1 - pi**-3 + 2 * pi**-5 - pi**-7 + pi**-9 - pi**-12 - 2 * pi**-14
    )
    E_C_neg = 2 * pi - pi**-1 + E_C_pos
    E[0] = pi**-12 + 2 * pi**-14
    if C > 0:
        E[0] = C * E_C_pos
    if C < 0:
        E[0] = -C * E_C_neg
    for l in GLUON_POWERS:
        E[2] += g2[l] * (2 * pi) ** l
    for n in FERMION_POWERS:
        E[1] -= g1[n] * (2 * pi) ** n
    for l in GLUON_POWERS:
        for n in FERMION_POWERS:
            if g2[l] != 0 and g1[n] != 0:
                E[3] += (l + n < 4) * (g2[l] > 0) * g2[l] * g1[n] * 2 * (2 * pi) ** (-l - n - 1)
                E[4] += (l + n < 4) * (g2[l] < 0) * g2[l] * g1[n] * 2 * (2 * pi) ** (-l - n)
                E[5] -= (l + n > 3) * g2[l] * g1[n] * 2 * (2 * pi) ** (-l - n - 1)
                E[6] += abs(g2[l] * g1[n]) * 2 * (2 * pi) ** -8
                g2[l] = 0
                g1[n] = 0
                break
            if g2[l] == 0 and g1[n] == 0:
                E[7] -= (2 * pi) ** (-l - n - 1)
                E[7] -= (2 * pi) ** (-l - n)
                break
    return E[0] + sum(E[1:8])


def paper_energies(coeffs: np.ndarray, base: float = 2 * math.pi) -> np.ndarray:
    """
    Vectorised paper formula. `coeffs` has shape (N, 7): i4,i3,i2,i1,i0,i-1,C
    (already halved, i.e. the actual half-integer values).

    `base` replaces 2π everywhere, including the π^-k constants (π := base/2),
    so a different base is a complete, self-consistent alternative theory.
    """
    pi = base / 2
    b = {k: base**k for k in range(-8, 5)}
    g = [coeffs[:, j].astype(float) for j in range(3)]
    f = [coeffs[:, 3 + j].astype(float) for j in range(3)]
    C = coeffs[:, 6].astype(float)

    E_C_pos = -pi + 2 * pi**-1 - pi**-3 + 2 * pi**-5 - pi**-7 + pi**-9 - pi**-12 - 2 * pi**-14
    E_C_neg = 2 * pi - pi**-1 + E_C_pos
    E0 = np.where(C > 0, C * E_C_pos, np.where(C < 0, -C * E_C_neg, pi**-12 + 2 * pi**-14))

    E2 = g[0] * b[4] + g[1] * b[3] + g[2] * b[2]
    E1 = -(f[0] * b[1] + f[1] * b[0] + f[2] * b[-1])
    rest = np.zeros(len(C))
    f_alive = [fv != 0 for fv in f]

    for gi, l in enumerate(GLUON_POWERS):
        gv = g[gi]
        g_nonzero = gv != 0
        active = np.ones(len(C), dtype=bool)
        for fi, n in enumerate(FERMION_POWERS):
            fv = f[fi]
            pair = active & g_nonzero & f_alive[fi]
            both_zero = active & ~g_nonzero & ~f_alive[fi]
            gf = gv * fv
            if l + n < 4:
                term = np.where(gv > 0, gf * 2 * b[-l - n - 1], gf * 2 * b[-l - n])
            else:
                term = -gf * 2 * b[-l - n - 1]
            term = term + np.abs(gf) * 2 * b[-8]
            rest = rest + np.where(pair, term, 0.0)
            rest = rest - np.where(both_zero, b[-l - n - 1] + b[-l - n], 0.0)
            f_alive[fi] = f_alive[fi] & ~pair
            active = active & ~(pair | both_zero)

    return E0 + E1 + E2 + rest


def coefficient_grid(j_gluon: int = 2, j_fermion: int = 3) -> np.ndarray:
    """All half-integer coefficient tuples: gluons in [-j_gluon, j_gluon], fermions in [-j_fermion, j_fermion], C in [-1, 1]."""
    gl = np.arange(-2 * j_gluon, 2 * j_gluon + 1) / 2
    fe = np.arange(-2 * j_fermion, 2 * j_fermion + 1) / 2
    cc = np.arange(-2, 3) / 2
    axes = np.meshgrid(gl, gl, gl, fe, fe, fe, cc, indexing="ij")
    return np.stack([a.ravel() for a in axes], axis=1)
