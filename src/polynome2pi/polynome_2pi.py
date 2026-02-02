import os
import argparse
import cmath
import datetime
from dataclasses import dataclass
from enum import Enum
from pathlib import Path

import numpy as np

import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection


RESULTS_DIR = Path("results")
RESULTS_DIR.mkdir(exist_ok=True)


class ScanSector(str, Enum):
    minimal = "minimal"
    light = "light"
    broad = "broad"
    nucleon = "nucleon"
    heavy = "heavy"


@dataclass
class PolynomeConfig:
    J4: int
    J3: int
    J2: int
    add_info: str = ""

    @property
    def name(self) -> str:
        suffix = f"_{self.add_info}" if self.add_info else ""
        return f"Polynom_{self.J4}{self.J3}{self.J2}{suffix}"


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description="Run P(2π) polynomial scan and plot results")

    parser.add_argument(
        "--sector",
        type=ScanSector,
        choices=list(ScanSector),
        required=True,
        help=(
            "Physical scan sector (polynomial depth). "
            f"Must be one of: {', '.join([e.value for e in ScanSector])}."
        ),
    )

    parser.add_argument(
        "--no-show",
        action="store_true",
        help="Do not open a GUI window; still save the PNG output",
    )

    parser.add_argument(
        "--debug",
        action="store_true",
        help="Print extra debug output (slow).",
    )

    return parser.parse_args(argv)


def select_preset_by_sector(sector: ScanSector) -> PolynomeConfig:
    if sector is ScanSector.minimal:
        return PolynomeConfig(J4=0, J3=0, J2=1, add_info=sector.value)

    if sector is ScanSector.light:
        return PolynomeConfig(J4=0, J3=1, J2=1, add_info=sector.value)

    if sector is ScanSector.broad:
        return PolynomeConfig(J4=1, J3=1, J2=2, add_info=sector.value)

    if sector is ScanSector.nucleon:
        return PolynomeConfig(J4=1, J3=1, J2=1, add_info=sector.value)

    if sector is ScanSector.heavy:
        return PolynomeConfig(J4=2, J3=2, J2=2, add_info=sector.value)

    raise ValueError(f"Unhandled sector: {sector}")


class PolynomeEngine:
    """
    Owns all precomputed constants and all scratch state used by energie().
    This eliminates globals while preserving performance (scratch reused).
    """

    def __init__(self, pow_off: int = 40):
        self.pi = cmath.pi
        self.two_pi = 2 * self.pi
        self.pow_off = pow_off

        # Precompute powers (2π)^k into a list for fast indexing
        self.two_pi_pow = [0.0] * (2 * pow_off + 1)
        for k in range(-pow_off, pow_off + 1):
            self.two_pi_pow[k + pow_off] = self.two_pi ** k

        # C-dependent constants
        pi = self.pi
        self.E_C_POS = (
            -pi
            + 2 * pi ** (-1)
            - pi ** (-3)
            + 2 * pi ** (-5)
            - pi ** (-7)
            + pi ** (-9)
            - pi ** (-12)
            - 2 * pi ** (-14)
        )
        self.E_C_NEG = 2 * pi - pi ** (-1) + self.E_C_POS
        self.E_C_ZERO = pi ** (-12) + 2 * pi ** (-14)

        # Scratch state (reused; no per-call allocations)
        self.E = [0.0] * 10
        self.g = [[0.0] * 10 for _ in range(10)]  # no aliasing

        # Hot precomputations for energie()
        off = self.pow_off
        p = self.two_pi_pow

        # (2π)^l for l in {4,3,2}
        self.POW_L_4 = p[4 + off]
        self.POW_L_3 = p[3 + off]
        self.POW_L_2 = p[2 + off]

        # (2π)^n for n in {1,0,-1}
        self.POW_N_1 = p[1 + off]
        self.POW_N_0 = p[0 + off]
        self.POW_N_M1 = p[-1 + off]

        # 2*(2π)^(-8) constant used in E6
        self.POW_NEG8_2 = 2.0 * p[-8 + off]

        # Precompute 9 combos for -l-n-1 and -l-n
        LS = (4, 3, 2)
        NS = (1, 0, -1)

        # include factor 2 for E3/E4/E5 terms
        self.POW_LN_M1 = [[0.0] * 3 for _ in range(3)]      # 2*(2π)^(-l-n-1)
        self.POW_LN_0 = [[0.0] * 3 for _ in range(3)]       # 2*(2π)^(-l-n)

        # exclude factor 2 for E7
        self.POW_LN_M1_ONLY = [[0.0] * 3 for _ in range(3)] # (2π)^(-l-n-1)
        self.POW_LN_0_ONLY = [[0.0] * 3 for _ in range(3)]  # (2π)^(-l-n)

        for li, l in enumerate(LS):
            for ni, n in enumerate(NS):
                self.POW_LN_M1[li][ni] = 2.0 * p[-l - n - 1 + off]
                self.POW_LN_0[li][ni] = 2.0 * p[-l - n + off]
                self.POW_LN_M1_ONLY[li][ni] = p[-l - n - 1 + off]
                self.POW_LN_0_ONLY[li][ni] = p[-l - n + off]

    def energie(self, i4, i3, i2, i1, i0, i_1, C):
        """
        Behavior-identical to the original Energie() nested-loop logic,
        but faster due to:
          - indexed pow tables (no dict)
          - local accumulators, fewer list writes
          - precomputed constants for the 9 (l,n) combos
        """
        g2 = self.g[2]
        g1 = self.g[1]

        # set inputs
        g2[4] = i4
        g2[3] = i3
        g2[2] = i2
        g1[1] = i1
        g1[0] = i0
        g1[-1] = i_1

        # C-dependent base term
        if C > 0:
            E0 = C * self.E_C_POS
        elif C < 0:
            E0 = -C * self.E_C_NEG
        else:
            E0 = self.E_C_ZERO

        # gluons and fermions
        E2 = g2[4] * self.POW_L_4 + g2[3] * self.POW_L_3 + g2[2] * self.POW_L_2
        E1 = -(g1[1] * self.POW_N_1 + g1[0] * self.POW_N_0 + g1[-1] * self.POW_N_M1)

        E3 = 0.0
        E4 = 0.0
        E5 = 0.0
        E6 = 0.0
        E7 = 0.0

        # exact original loop/break semantics
        for li, l in enumerate((4, 3, 2)):
            for ni, n in enumerate((1, 0, -1)):
                gl = g2[l]
                gn = g1[n]

                if gl != 0 and gn != 0:
                    ln = l + n

                    if ln < 4:
                        if gl > 0:
                            E3 += gl * gn * self.POW_LN_M1[li][ni]
                        else:
                            E4 += gl * gn * self.POW_LN_0[li][ni]

                    if ln > 3:
                        E5 -= gl * gn * self.POW_LN_M1[li][ni]

                    prod = gl * gn
                    E6 += (prod if prod >= 0 else -prod) * self.POW_NEG8_2

                    g2[l] = 0
                    g1[n] = 0
                    break

                if gl == 0 and gn == 0:
                    E7 -= self.POW_LN_M1_ONLY[li][ni]
                    E7 -= self.POW_LN_0_ONLY[li][ni]
                    break

        total = E0 + E1 + E2 + E3 + E4 + E5 + E6 + E7

        # preserve side-effects (main loop reads E[0])
        E = self.E
        E[0] = total
        E[1] = E1
        E[2] = E2
        E[3] = E3
        E[4] = E4
        E[5] = E5
        E[6] = E6
        E[7] = E7

        return total


# ------------------------- extracted helpers -------------------------

def build_colors():
    # colors (index 1..26 used)
    return [
        "#FFFFFF", "#000000", "#F60000", "#05FB4F", "#CFCF00", "#000000", "#07FCE4",
        "#F700D2", "#00F73E", "#7BB91F", "#A9BF06", "#047619", "#047619", "#789E20",
        "#CFCF00", "#CF00B7", "#EC61A9", "#FA9805", "#4200F6", "#495999", "#B91F50",
        "#CB4088", "#F90404", "#000000", "#4D8E2F", "#499999", "#F50606", "#146108",
    ]


def build_particle_table():
    Obj = [
        ["Name", "m_e", "E", "-SD", "+SD", "Halbwertszeit T in sec", "Charge", "Spin", "P.", "name"],
        ["e", "1.00000000000(31)", "1.00000000000", "-0.005", "0.000", " ", "-1", "1/2", "", "e"],
        ["u", "4.18(-0.51)(0.96)", "4.18", "-0.51", "0.96", " ", "+2/3", "", "", "u"],
        ["d", "9.14(-0.33)(0.94)", "9.14", "-0.33", "0.94", " ", "-1/3", "", "", "d"],
        ["s", "182.8(-6.6)(16.8)", "182.8", "-6.6", "16.8", " ", "-1/3", "", "", "s"],
        ["Muon", "206.7682827(46)", "206.7682827", "-0.0000046", "0.0000046", "2.1969811(22) e -6", "0", "1", "", "muon"],
        ["Pion 0", "264.1430(9)", "264.1430", "-0.0009", "0.0009", "8.52(18) e -17", "0", "0", "-", r"$u\overline{d}-\overline{u}d$"],
        ["Pion +-", "273.13243(35)", "273.13243", "-0.00035", "0.00035", "2.6033(5) e -8", "+-1", "0", "-", r"$u\overline{u},\overline{d}d$"],
        ["K +-", "966.102(21)", "966.102", "-0.021", "0.021", "1.2380(20) e -8", "+-1", "0", "-", r"$u\overline{s},s\overline{u}$"],
        ["KL 0", "973.800(26)", "973.800", "-0.026", "0.026", "5.116(21) e -8", "0", "0", "-", r"$d\overline{s},s\overline{d}$  "],
        ["KS 0", "973.800(26)", "973.800", "-0.026", "0.026", "8.954(4) e -11", "0", "", "", r"$d\overline{s},s\overline{d}$"],
        ["Eta", "1072.139(35)", "1072.139", "-0.035", "0.035", "5 e -19", "0", "0", "-", r"$u\overline{u}+\overline{d}d-2s\overline{s}$"],
        ["Rho +-", "1506(1)", "1506", "-1", "1", "4 e -24", "-+1", "1", "-", r"$u\overline{u},\overline{d}d$"],
        ["Rho 0", "1517.14(49)", "1517.14", "-0.49", "0.49", "4 e -24", "0", "1", "-", r"$u\overline{u}-\overline{d}d$"],
        ["Omega", "1531.62(25)", "1531.62", "-0.25", "0.25", "7.75(7) e -23", "0", "1", "-", r"$u\overline{u}+\overline{d}d$"],
        ["K* +-", "1745.2(1)", "1745.2", "-0.1", "0.1", "1.3 e -23", "+-1", "", "", r"$d\overline{s},s\overline{d}$"],
        ["K* 0", "1752.6(1)", "1752.6", "-0.1", "0.1", "1.3 e -23", "0", "", "", r"$d\overline{s},s\overline{d}$"],
        ["Proton", "1836.152673426(32)", "1836.152673426", "-0.000000032", "0.000000032", " ", "1", "1/2", "1", "uud"],
        ["H", "1837.47(-0.29)(0.20)", "1837.47", "-0.29", "0.20", " ", "0", "", "", "H"],
        ["Neutron", "1838.68366200(74)", "1838.68366200", "-0.00000074", "0.00000074", "878.4(5)", "0", "1/2", "1", "udd"],
        ["Eta`", "1874.32(11)", "1874.32", "-0.11", "0.11", "3.32(15) e -21", "0", "0", "-", r"$u\overline{u}+\overline{d}d+s\overline{s}$"],
        ["Phi", "1995.035(31)", "1995.035", "-0.031", "0.031", "1.55(0,01) e -22", "0", "1", "-", r"$s\overline{s}(most)$"],
        ["c", "2485(-39)(39)", "2485", "-39", "39", " ", "+2/3", "", "", "c", 0],
        ["Tau", "3477.23(23)", "3477.23", "-0.23", "0.23", "290.3(5) e -15", "-1", "1/2", "", "tau"],
        ["D 0", "3649.38(10)", "3649.38", "-0.10", "0.10", "4.101(15) e -13", "+-1", "0", "-", r"$c\overline{u},u\overline{c}$"],
        ["D +", "3658.81(10)", "3658.81", "-0.10", "0.10", "1.040(7) e -12", "+-1", "0", "-", r"$c\overline{d},d\overline{c}$"],
        ["Deuteron", "3670.4829677(11)", "3670.4829677", "-0.0000011", "0.0000011", " ", "0", "", "", "D", -5.5],
        ["DS +", "3851.94(13)", "3851.94", "-0.13", "0.13", "5.04(4) e -13", "+-1", "0", "-", r"$c\overline{s},s\overline{c}$"],
        ["Higgs", "244830(210)", "244830", "-210", "210", " ", "0", "0", "", "Higgs"],
        ["t", "337710(570)", "337710", "-570", "570", " ", "+2/3", "", "", "t"],
    ]

    # numeric helpers for matching (avoid float conversion in hot loop)
    obj_E = [0.0] * 27
    obj_min = [0.0] * 27
    obj_max = [0.0] * 27
    for jj in range(1, 27):
        obj_E[jj] = float(Obj[jj][2])
        obj_min[jj] = float(Obj[jj][3])
        obj_max[jj] = float(Obj[jj][4])

    return Obj, obj_E, obj_min, obj_max


def allocate_result_arrays():
    D_i_N = [0.0] * 300
    D_i_c_ = [""] * 35
    Cnt = [0] * 520

    Emax = np.zeros((35, 520), dtype=float)
    Emin = np.zeros((35, 520), dtype=float)
    i_Emax = np.zeros((35, 520), dtype=int)
    i_Emin = np.zeros((35, 520), dtype=int)
    Dmax = np.zeros((35, 520, 7), dtype=int)
    Dmin = np.zeros((35, 520, 7), dtype=int)

    return D_i_N, D_i_c_, Cnt, Emax, Emin, i_Emax, i_Emin, Dmax, Dmin


def init_plot_buffers():
    xs_by_j = [[] for _ in range(27)]
    ys_by_j = [[] for _ in range(27)]
    grey_segments = []
    return xs_by_j, ys_by_j, grey_segments


def run_scan(
    *,
    engine: PolynomeEngine,
    config: PolynomeConfig,
    is_heavy: bool,
    is_nucleon: bool,
    debug: bool,
    obj_E,
    obj_min,
    obj_max,
    Cnt,
    Emax,
    Emin,
    i_Emax,
    i_Emin,
    Dmax,
    Dmin,
    xs_by_j,
    ys_by_j,
    grey_segments,
):
    """
    Executes the nested scan loops. Returns (i_T, i_T1, mmax).
    Mutates the provided arrays/buffers.
    """
    energie = engine.energie
    E_local = engine.E

    i_T = 0
    i_T1 = 0
    mmax = 0
    ct = 0

    for i4 in range(-2 * config.J4, 2 * config.J4 + 1):
        i4h = 0.5 * i4
        for i3 in range(-2 * config.J3, 2 * config.J3 + 1):
            i3h = 0.5 * i3
            for i2 in range(-2 * config.J2, 2 * config.J2 + 1):
                i2h = 0.5 * i2

                if debug:
                    print("i4", i4, "i3", i3, "i2", i2, "i_T1", i_T1)

                for i1 in range(-6, 7):
                    i1h = 0.5 * i1
                    for i0 in range(-6, 7):
                        i0h = 0.5 * i0
                        for i_1 in range(-6, 7):
                            i_1h = 0.5 * i_1
                            for C in range(-2, 3):
                                Ch = 0.5 * C

                                energie(i4h, i3h, i2h, i1h, i0h, i_1h, Ch)
                                E0 = E_local[0]
                                if E0 < 0:
                                    continue

                                m = int(256 + 32 * i4 + 4 * i3 + i2)
                                if m > mmax:
                                    mmax = m
                                    ct = 0
                                ct += 1
                                Cnt[mmax] = ct

                                if is_heavy and E0 < 1500:
                                    continue
                                if is_nucleon and (E0 < 1836 or E0 > 1839):
                                    continue

                                i_T += 1
                                flag_match = 0

                                for j in range(1, 27):
                                    Ej = obj_E[j]
                                    if (E0 - Ej <= obj_max[j]) and (E0 - Ej >= obj_min[j]):
                                        i_T1 += 1

                                        if Emax[j, m] <= E0:
                                            Emax[j, m] = E0
                                            i_Emax[j, m] = i_T
                                            Dmax[j, m, 0] = i4
                                            Dmax[j, m, 1] = i3
                                            Dmax[j, m, 2] = i2
                                            Dmax[j, m, 3] = i1
                                            Dmax[j, m, 4] = i0
                                            Dmax[j, m, 5] = i_1
                                            Dmax[j, m, 6] = C

                                        if Emin[j, m] >= E0 or Emin[j, m] == 0:
                                            Emin[j, m] = E0
                                            i_Emin[j, m] = i_T
                                            Dmin[j, m, 0] = i4
                                            Dmin[j, m, 1] = i3
                                            Dmin[j, m, 2] = i2
                                            Dmin[j, m, 3] = i1
                                            Dmin[j, m, 4] = i0
                                            Dmin[j, m, 5] = i_1
                                            Dmin[j, m, 6] = C

                                        xs_by_j[j].append(i_T)
                                        ys_by_j[j].append(E0)
                                        flag_match = 1

                                if E0 > 0 and flag_match == 0:
                                    grey_segments.append([(i_T, E0), (i_T + 1, E0)])

    return i_T, i_T1, mmax


def draw_points(xs_by_j, ys_by_j, grey_segments, colors):
    for j in range(1, 27):
        if xs_by_j[j]:
            plt.scatter(xs_by_j[j], ys_by_j[j], s=80, c=colors[j], marker=".", linewidths=0)

    if grey_segments:
        lc = LineCollection(grey_segments, colors="#C0BCBC", linewidths=1)
        plt.gca().add_collection(lc)


def label_offsets_for_sector(sector: ScanSector, j: int):
    """
    Returns (X, Y, fs) where X is the whole list so we can use X[j].
    Keeping this identical to your working version.
    """
    if sector is ScanSector.broad:
        X = [0, 4, 7, 9, 6, 9, 5, 13, -21, -12, 5, -19, 6, -28, -17, -21, -13, -20, -12, -10, -37, -9, -5.5, -22, -11, 0, 0, 0]
        return X, -20, 12

    if sector is ScanSector.light:
        X = [0, 1.2, 2.5, 3, 1, 1, 0.7, 1]
        return X, -4, 16

    if sector is ScanSector.minimal:
        X = [0, 0.2, 0.2, 0.2]
        return X, -1, 16

    if sector is ScanSector.nucleon:
        X = [0, 2, 4, 5, 3, 4, 4, 8, -13, -7, 2, -11, 2, -14, -7, -11, -7, -12, -8, -3, -0, -0, -0, -22, -11, 0, 0, 0]
        return X, -0, 16

    if sector is ScanSector.heavy:
        X = [0, 5, 10, 15, 20, 5, 9, 30, -35, -20, 7, -40, 25, 15, 4, 35, 20, -13, -10, -7, 20, 20, 10, -10, -10, -20, 0, 0, 0]
        X[j] *= 2
        return X, -50, 12

    raise ValueError(f"Unhandled sector: {sector}")


def write_report_and_labels(
    *,
    f,
    sector: ScanSector,
    Obj,
    colors,
    D_i_N,
    D_i_c_,
    Cnt,
    Emax,
    Emin,
    i_Emax,
    i_Emin,
    Dmax,
    Dmin,
):
    print("..............   plotting   .................")
    # (i_T / i_T1 prints are done in main because it has those values)

    print(
        "{0:10}{1:5}{2:20}{3:16}{4:8}{5:5}{6:5}{7:5}{8:5}{9:5}{10:5}{11:5}".format(
            "particle", "", "", "     theory: E", "   total ", "  i4", "  i3", "  i2", "  i1", "  i0", " i-1", "  C"
        ),
        file=f,
    )

    for j in range(1, 29):
        m_min = float(Obj[j][2]) + float(Obj[j][3])
        m_max = float(Obj[j][2]) + float(Obj[j][4])
        p_g = len(Obj[j][2])
        m_min = float(str(m_min)[:p_g])
        m_max = float(str(m_max)[:p_g])
        flag = 0

        for m in range(1, 512):
            if i_Emax[j, m] == 0:
                continue
            if flag == 0:
                print(Obj[j][0], file=f)
                flag += 1

            E_mean = (Emax[j, m] + Emin[j, m]) / 2
            Di_E = (i_Emax[j, m] - i_Emin[j, m] + 1)
            i_Emax[j, 0] += abs(Di_E)
            Cnt_ = Cnt[m]
            D_i_c = round((abs(Di_E)) * 100 / Cnt_, 5)
            i_Emax[j, 516] += Cnt_
            D_i_N[j] = float(i_Emax[j, 0]) * 100 / Cnt[m]

            Emax[j, m] = float(str(Emax[j, m])[:p_g])
            Emin[j, m] = float(str(Emin[j, m])[:p_g])
            E_mean = float(str(E_mean)[:p_g])

            print(
                "{0:10}{1:5}{2:20}{3:16}{4:8}{5:5}{6:5}{7:5}{8:5}{9:5}{10:5}{11:5}".format(
                    "", "max   ", m_max, Emax[j, m], i_Emax[j, m],
                    Dmax[j, m, 0] / 2, Dmax[j, m, 1] / 2, Dmax[j, m, 2] / 2,
                    Dmax[j, m, 3] / 2, Dmax[j, m, 4] / 2, Dmax[j, m, 5] / 2,
                    Dmax[j, m, 6] / 2,
                ),
                file=f,
            )
            print(
                "{0:10}{1:5}{2:20}{3:16}{4:8}{5:9}{6:7}{7:2}{8:12}{9:7}{10:2}".format(
                    "", "mean  ", Obj[j][1], E_mean, Di_E, "", "", "", "", "", ""
                ),
                file=f,
            )
            print(
                "{0:10}{1:5}{2:20}{3:16}{4:8}{5:5}{6:5}{7:5}{8:5}{9:5}{10:5}{11:5}".format(
                    "", "min   ", m_min, Emin[j, m], i_Emin[j, m],
                    Dmin[j, m, 0] / 2, Dmin[j, m, 1] / 2, Dmin[j, m, 2] / 2,
                    Dmin[j, m, 3] / 2, Dmin[j, m, 4] / 2, Dmin[j, m, 5] / 2,
                    Dmin[j, m, 6] / 2,
                ),
                file=f,
            )
            print(
                "{0:10}{1:5}{2:20}{3:15}{4:8}{5:5}{6:8}{7:9}{8:10}{9:2}".format(
                    "", "", "", "         ∆ abs(i)", abs(Di_E),
                    "  Cts", Cnt_, "  ∆i/(2pi)", D_i_c, " %"
                ),
                file=f,
            )

            # Restore labels in the plot
            if flag == 1:
                X, Y, fs = label_offsets_for_sector(sector, j)
                plt.text(
                    i_Emin[j, m] + 10000 * X[j],
                    E_mean + Y,
                    Obj[j][9],
                    fontsize=fs,
                    color=colors[j],
                )
                flag = 2

        if j > 1 and m == 511 and i_Emax[j, 516] > 0:
            D_i_c = round(float(i_Emax[j, 0]) * 100 / i_Emax[j, 516], 5)
            D = str(D_i_c)
            D_i_c_[j] = "".rjust(10 - len(D), " ") + D + " %"
            print(
                "{0:10}{1:5}{2:20}{3:15}{4:8}{5:5}{6:8}{7:9}{8:10}{9:2}".format(
                    "", "total", "", "             Σ ∆i", i_Emax[j, 0],
                    "  Cts", i_Emax[j, 516], "  ∆i/(2pi)", D_i_c, " %"
                ),
                file=f,
            )

        if flag == 0:
            print(
                "{0:10}{1:5}{2:20}{3:16}{4:8}".format(
                    Obj[j][0], "mean  ", Obj[j][1], "  ", " only with i4 > 1"
                ),
                file=f,
            )


def add_reference_lines(sector: ScanSector, i_T: int):
    x_a = 0
    x_m = i_T * 1 / 5

    plt.ylabel("Energy in $m_e$")
    plt.xlabel("N")

    TWO_PI = float(2 * cmath.pi)
    i2_ = TWO_PI ** 2
    i3_ = TWO_PI ** 3
    i4_ = TWO_PI ** 4
    i5_ = 0.5 * (i4_ + i3_ + i2_)
    i6_ = i4_ + i3_ + i2_
    i7_ = 2.5 * i4_ - 1.5 * i3_ - 0.5 * i2_
    i8_ = 1.5 * i4_ + i3_ + i2_
    i10_ = 1.5 * i4_ + 0.5 * i3_ - 0.5 * i2_

    if sector is ScanSector.broad:
        plt.plot([x_a, x_m], [i4_, i4_], "k", linewidth=1)
        plt.text(x_a, i4_ + 15, r"$(2\pi)^4$", fontsize=12, color="blue")
        plt.plot([x_a, x_m], [i3_, i3_], "k", linewidth=1)
        plt.text(x_a, i3_ + 15, r"$(2\pi)^3$", fontsize=12, color="blue")
        plt.plot([x_a, x_m], [i2_, i2_], "k", linewidth=1)
        plt.text(x_a, i2_ + 15, r"$(2\pi)^2$", fontsize=12, color="blue")
        plt.plot([x_a, x_m], [i5_, i5_], "k", linewidth=1)
        plt.text(x_a, i5_ + 15, r"$1/2((2\pi)^4+(2\pi)^3+(2\pi)^2)$", fontsize=12, color="blue")
        plt.plot([x_a, x_m], [i6_, i6_], "k", linewidth=1)
        plt.text(x_a, i6_ + 15, r"$(2\pi)^4+(2\pi)^3+(2\pi)^2$", fontsize=12, color="blue")

    if sector is ScanSector.heavy:
        plt.plot([x_a, x_m], [i7_, i7_], "k", linewidth=1)
        plt.text(x_a, i7_ + 15, r"$5/2(2\pi)^4-3/2(2\pi)^3-1/2(2\pi)^2$", fontsize=12, color="blue")
        plt.plot([x_a, 2 * x_m], [i8_, i8_], "k", linewidth=1)
        plt.text(x_a, i8_ + 15, r"$3/2(2\pi)^4+(2\pi)^3+(2\pi)^2$", fontsize=12, color="blue")
        plt.plot([x_a, 3 * x_m], [i10_, i10_], "k", linewidth=1)
        plt.text(x_a, i10_ + 15, r"$3/2(2\pi)^4+1/2(2\pi)^3-1/2(2\pi)^2$", fontsize=12, color="blue")

    if sector in (ScanSector.broad, ScanSector.light):
        plt.plot([x_a, x_m], [i3_, i3_], "k", linewidth=1)
        plt.text(x_a, i3_ + 15, r"$(2\pi)^3$", fontsize=12, color="blue")

    if sector in (ScanSector.broad, ScanSector.light, ScanSector.minimal):
        plt.plot([x_a, x_m], [i2_, i2_], "k", linewidth=1)
        plt.text(x_a, i2_ + 15, r"$(2\pi)^2$", fontsize=12, color="blue")

    # x-limits similar to your previous plots
    if sector is ScanSector.minimal:
        plt.xlim(-1000, i_T + 10000)
    elif sector is ScanSector.nucleon:
        plt.xlim(0, i_T)
    else:
        plt.xlim(-10000, i_T + 30000)


def add_legend_panels(sector: ScanSector, i_T: int, Obj, colors, i_Emax, D_i_c_):
    if sector is ScanSector.broad:
        x_a = i_T * 0.65
        dx = i_T * 0.08
        i = -20
        for jj in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 16, 17, 18, 19]:
            particle = str(Obj[jj][9])
            plt.text(x_a, i, Obj[jj][0])
            plt.text(x_a + dx, i, particle, color=colors[jj])
            plt.text(x_a + 2 * dx, i, i_Emax[jj, 0])
            plt.text(x_a + 3 * dx, i, D_i_c_[jj])
            i += 100
        plt.text(x_a + 2 * dx, i, "  ∆i  ")
        plt.text(x_a + 3 * dx, i, " ∆i/(2pi) ")

    if sector is ScanSector.light:
        x_a = i_T * 0.70
        dx = i_T * 0.10
        i = 0
        for jj in [1, 2, 3, 4, 5, 6, 7]:
            particle = str(Obj[jj][9])
            plt.text(x_a, i, Obj[jj][0])
            plt.text(x_a + dx, i, particle, color=colors[jj])
            plt.text(x_a + 2 * dx, i, i_Emax[jj, 0])
            plt.text(x_a + 3 * dx, i, D_i_c_[jj])
            i += 20
        plt.text(x_a + 2 * dx, i, "  ∆i  ")
        plt.text(x_a + 3 * dx, i, " ∆i/(2pi) ")

    if sector is ScanSector.minimal:
        x_a = i_T * 0.89
        dx = i_T * 0.10
        i = 0
        for jj in [1, 2, 3, 4, 5, 6, 7]:
            particle = str(Obj[jj][9])
            plt.text(x_a, i, Obj[jj][0])
            plt.text(x_a + dx, i, particle, color=colors[jj])
            plt.text(x_a + 2 * dx, i, i_Emax[jj, 0])
            plt.text(x_a + 3 * dx, i, D_i_c_[jj])
            i += 5
        plt.text(x_a + 2 * dx, i, "  ∆i  ")
        plt.text(x_a + 3 * dx, i, " ∆i/(2pi) ")

    if sector is ScanSector.nucleon:
        fs = 14
        x_a = i_T * 0.70
        dx = i_T * 0.10
        i = 0
        m_H = 1836.152673426 + 1
        plt.plot([1, i_T], [m_H, m_H], "k", linewidth=1)
        plt.text(1, 1837, "$m_{Proton} + m_e$", fontsize=fs, color="blue")
        for jj in [17, 19]:
            plt.text(1050, float(Obj[jj][2]), Obj[jj][0], fontsize=fs, color=colors[jj])
            particle = str(Obj[jj][9])
            plt.text(x_a, i, Obj[jj][0])
            plt.text(x_a + dx, i, particle, color=colors[jj])
            plt.text(x_a + 2 * dx, i, i_Emax[jj, 0])
            plt.text(x_a + 3 * dx, i, D_i_c_[jj])
            i += 20

    if sector is ScanSector.heavy:
        x_a = i_T * 0.65
        dx = i_T * 0.085
        i = 1500
        for jj in [15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25]:
            particle = str(Obj[jj][9])
            plt.text(x_a, i, Obj[jj][0])
            plt.text(x_a + dx, i, particle, color=colors[jj])
            plt.text(x_a + 2 * dx, i, i_Emax[jj, 0])
            plt.text(x_a + 3 * dx, i, D_i_c_[jj])
            i += 100
        plt.text(x_a + 2 * dx, i, "  ∆i  ")
        plt.text(x_a + 3 * dx, i, " ∆i/(2pi) ")


# ------------------------- main orchestration -------------------------

def main(argv=None):
    start_time_stamp = datetime.datetime.now()

    args = parse_args(argv)
    sector = args.sector
    config = select_preset_by_sector(sector)
    no_show = args.no_show
    DEBUG = args.debug

    # hoist sector checks (avoid Enum comparisons in the hot loop)
    is_heavy = (sector is ScanSector.heavy)
    is_nucleon = (sector is ScanSector.nucleon)

    engine = PolynomeEngine()

    # data
    colors = build_colors()
    Obj, obj_E, obj_min, obj_max = build_particle_table()
    D_i_N, D_i_c_, Cnt, Emax, Emin, i_Emax, i_Emin, Dmax, Dmin = allocate_result_arrays()
    xs_by_j, ys_by_j, grey_segments = init_plot_buffers()

    # output paths
    txt_path = RESULTS_DIR / f"{config.name}.txt"
    png_path = RESULTS_DIR / f"{config.name}.png"

    with open(txt_path, "w", encoding="utf8") as f:
        # scan
        i_T, i_T1, _mmax = run_scan(
            engine=engine,
            config=config,
            is_heavy=is_heavy,
            is_nucleon=is_nucleon,
            debug=DEBUG,
            obj_E=obj_E,
            obj_min=obj_min,
            obj_max=obj_max,
            Cnt=Cnt,
            Emax=Emax,
            Emin=Emin,
            i_Emax=i_Emax,
            i_Emin=i_Emin,
            Dmax=Dmax,
            Dmin=Dmin,
            xs_by_j=xs_by_j,
            ys_by_j=ys_by_j,
            grey_segments=grey_segments,
        )

        # plot points (batched)
        draw_points(xs_by_j, ys_by_j, grey_segments, colors)

        # report + particle labels
        print("possible ET:", i_T, "real ET:", i_T1)
        write_report_and_labels(
            f=f,
            sector=sector,
            Obj=Obj,
            colors=colors,
            D_i_N=D_i_N,
            D_i_c_=D_i_c_,
            Cnt=Cnt,
            Emax=Emax,
            Emin=Emin,
            i_Emax=i_Emax,
            i_Emin=i_Emin,
            Dmax=Dmax,
            Dmin=Dmin,
        )

        # reference lines and legend panels
        add_reference_lines(sector, i_T)
        add_legend_panels(sector, i_T, Obj, colors, i_Emax, D_i_c_)

    # save/show
    fig = plt.gcf()
    fig.set_size_inches(10, 6)
    fig.savefig(png_path, dpi=100)

    end_time_stamp = datetime.datetime.now()
    delta = end_time_stamp - start_time_stamp
    print(f"took {delta.seconds} seconds")
    print(f"wrote: {txt_path}")
    print(f"wrote: {png_path}")

    if not no_show:
        plt.show()


if __name__ == "__main__":
    main()