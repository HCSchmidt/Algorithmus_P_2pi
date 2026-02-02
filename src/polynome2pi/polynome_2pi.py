import os
import argparse
import cmath
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import matplotlib.colors as mcolors
import datetime

RESULTS_DIR = Path("results")
RESULTS_DIR.mkdir(exist_ok=True)

pi = cmath.pi

# ---- Precomputed constants for performance ----
TWO_PI = 2 * pi

# Precompute powers of (2*pi) used in Energie(). Keep a generous range.
TWO_PI_POW = {k: (TWO_PI ** k) for k in range(-30, 31)}

# Precompute the C-dependent energy constants (previously recomputed every call).
E_C_POS = (
    -pi
    + 2 * pi ** (-1)
    - pi ** (-3)
    + 2 * pi ** (-5)
    - pi ** (-7)
    + pi ** (-9)
    - pi ** (-12)
    - 2 * pi ** (-14)
)
E_C_NEG = 2 * pi - pi ** (-1) + E_C_POS
E_C_ZERO = pi ** (-12) + 2 * pi ** (-14)

E = [0] * 10

# FIX 1) avoid shared-list aliasing (each row independent)
g = [[0] * 10 for _ in range(10)]


@dataclass
class PolynomeConfig:
    J4: int
    J3: int
    J2: int
    add_info: str = ""

    @property
    def name(self) -> str:
        """Human-readable name used for output files."""
        suffix = f"_{self.add_info}" if self.add_info else ""
        return f"Polynom_{self.J4}{self.J3}{self.J2}{suffix}"


def Energie(i4, i3, i2, i1, i0, i_1, C):
    g[2][4] = i4
    g[2][3] = i3
    g[2][2] = i2
    g[1][1] = i1
    g[1][0] = i0
    g[1][-1] = i_1
    E[0] = 0
    E[1] = 0
    E[2] = 0
    E[3] = 0
    E[4] = 0
    E[5] = 0
    E[6] = 0
    E[7] = 0

    # C-dependent term (now precomputed constants)
    E[0] = E_C_ZERO
    if C > 0:
        E[0] = C * E_C_POS
    if C < 0:
        E[0] = -C * E_C_NEG

    for l in range(4, 1, -1):  # Gluonen r b g
        E[2] += g[2][l] * TWO_PI_POW[l]
    for n in range(1, -2, -1):  # e, u, d
        E[1] -= g[1][n] * TWO_PI_POW[n]

    for l in range(4, 1, -1):
        for n in range(1, -2, -1):
            if g[2][l] != 0 and g[1][n] != 0:
                E[3] += (l + n < 4) * (g[2][l] > 0) * g[2][l] * g[1][n] * 2 * TWO_PI_POW[-l - n - 1]
                E[4] += (l + n < 4) * (g[2][l] < 0) * g[2][l] * g[1][n] * 2 * TWO_PI_POW[-l - n]
                E[5] -= (l + n > 3) * g[2][l] * g[1][n] * 2 * TWO_PI_POW[-l - n - 1]
                E[6] += abs(g[2][l] * g[1][n]) * 2 * TWO_PI_POW[-8]
                g[2][l] = 0
                g[1][n] = 0
                break
            if g[2][l] == 0 and g[1][n] == 0:
                E[7] -= TWO_PI_POW[-l - n - 1]
                E[7] -= TWO_PI_POW[-l - n]
                break

    E[0] += E[1] + E[2] + E[3] + E[4] + E[5] + E[6] + E[7]
    return E[0]


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description="Run P(2π) polynomial scan and plot results")
    parser.add_argument("--op", type=int, default=2, help="Plot/scan preset option (1-5). Default: 2")
    parser.add_argument("--no-show", action="store_true", help="Do not open a GUI window; still save the PNG output")
    args = parser.parse_args(argv)

    op = args.op
    no_show = args.no_show

    return op, no_show


def select_preset(op: int):
    if op == 1:  # slow (20 minutes)
        config = PolynomeConfig(J4=1, J3=1, J2=2)
    elif op == 2:  # fast, for u , d , s
        config = PolynomeConfig(J4=0, J3=1, J2=1)
    elif op == 3:  # fast, for u , d
        config = PolynomeConfig(J4=0, J3=0, J2=1)
    elif op == 4:  # fast, for Proton H-Atom Neutron
        config = PolynomeConfig(J4=1, J3=1, J2=1, add_info="Atom")
    elif op == 5:  #  E<1500! slow, more particles with c,tau
        config = PolynomeConfig(J4=2, J3=2, J2=2, add_info="Tau")
    else:
        raise ValueError("op must be 1..5")

    return config


def main(argv=None):

    start_time_stamp = datetime.datetime.now()
    Op, no_show = parse_args(argv)
    config = select_preset(Op)

    DEBUG = False

    with open(os.path.join(RESULTS_DIR, f"{config.name}.txt"), "w", encoding="utf8") as f:
        i_T = 0
        i_T1 = 0

        i1 = 0
        i0 = 0
        i_1 = 0
        mmax = 0
        m = 0

        j = 0
        D_i_N = [0] * 300

        # FIX 1) avoid shared-list aliasing for big matrices too
        Obj = [[0] * 13 for _ in range(35)]
        N_E = [0] * 100000
        N_T = [[0] * 100 for _ in range(100)]
        E_t = [0] * 100000

        X = [0] * 35
        D_i_c_ = [0] * 35
        Cnt = [0] * 520
        Emax = np.zeros((35, 520), dtype=float)
        Emin = np.zeros((35, 520), dtype=float)
        i_Emax = np.zeros((35, 520), dtype=int)
        i_Emin = np.zeros((35, 520), dtype=int)
        Dmax = np.zeros((35, 520, 7), dtype=int)
        Dmin = np.zeros((35, 520, 7), dtype=int)

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

        # Precompute numeric values for fast matching (avoid float() in the inner loop)
        obj_E = [0.0] * 27
        obj_min = [0.0] * 27
        obj_max = [0.0] * 27
        for jj in range(1, 27):
            obj_E[jj] = float(Obj[jj][2])
            obj_min[jj] = float(Obj[jj][3])
            obj_max[jj] = float(Obj[jj][4])

        F = [
            "#FFFFFF",
            "#000000",
            "#F60000",
            "#05FB4F",
            "#CFCF00",
            "#000000",
            "#07FCE4",
            "#F700D2",
            "#00F73E",
            "#7BB91F",
            "#A9BF06",
            "#047619",
            "#047619",
            "#789E20",
            "#CFCF00",
            "#CF00B7",
            "#EC61A9",
            "#FA9805",
            "#4200F6",
            "#495999",
            "#B91F50",
            "#CB4088",
            "#F90404",
            "#000000",
            "#4D8E2F",
            "#499999",
            "#F50606",
            "#146108",
        ]

        # --- Plot batching for performance ---
        xs_by_j = {jj: [] for jj in range(1, 27)}
        ys_by_j = {jj: [] for jj in range(1, 27)}
        grey_segments = []

        # FIX 2) local bindings for faster inner-loop access (locals are faster than globals)
        E_local = E
        Emax_local = Emax
        Emin_local = Emin
        i_Emax_local = i_Emax
        i_Emin_local = i_Emin
        Dmax_local = Dmax
        Dmin_local = Dmin
        obj_E_local = obj_E
        obj_min_local = obj_min
        obj_max_local = obj_max
        xs_by_j_local = xs_by_j
        ys_by_j_local = ys_by_j
        grey_segments_local = grey_segments
        Energie_func = Energie  # local bind function lookup too

        for i5 in [0]:
            for i4 in range(-2 * config.J4, 2 * config.J4 + 1):
                for i3 in range(-2 * config.J3, 2 * config.J3 + 1):
                    for i2 in range(-2 * config.J2, 2 * config.J2 + 1):
                        if DEBUG:
                            print("i4", i4, "i3", i3, "i2", i2, "i1", i1, "i_T1", i_T1)

                        for i1 in range(-6, 7):
                            for i0 in range(-6, 7):
                                for i_1 in range(-6, 7):
                                    for C in range(-2, 3):
                                        Energie_func(i4 / 2, i3 / 2, i2 / 2, i1 / 2, i0 / 2, i_1 / 2, C / 2)

                                        # FIX 3) read E0 once per combination
                                        E0 = E_local[0]
                                        if E0 < 0:
                                            continue

                                        m = int(256 + 32 * i4 + 4 * i3 + i2)
                                        if m > mmax:
                                            mmax = m
                                            ct = 0
                                        ct += 1
                                        Cnt[mmax] = ct

                                        if Op == 5 and E0 < 1500:
                                            continue
                                        if Op == 4 and (E0 < 1836 or E0 > 1839):
                                            continue

                                        i_T += 1
                                        flag = 0

                                        for j in range(1, 27):
                                            Ej = obj_E_local[j]
                                            if (E0 - Ej <= obj_max_local[j]) and (E0 - Ej >= obj_min_local[j]):
                                                i_T1 += 1
                                                if Emax_local[j, m] <= E0:
                                                    Emax_local[j, m] = E0
                                                    i_Emax_local[j, m] = i_T
                                                    Dmax_local[j, m, 0] = i4
                                                    Dmax_local[j, m, 1] = i3
                                                    Dmax_local[j, m, 2] = i2
                                                    Dmax_local[j, m, 3] = i1
                                                    Dmax_local[j, m, 4] = i0
                                                    Dmax_local[j, m, 5] = i_1
                                                    Dmax_local[j, m, 6] = C
                                                if Emin_local[j, m] >= E0 or Emin_local[j, m] == 0:
                                                    Emin_local[j, m] = E0
                                                    i_Emin_local[j, m] = i_T
                                                    Dmin_local[j, m, 0] = i4
                                                    Dmin_local[j, m, 1] = i3
                                                    Dmin_local[j, m, 2] = i2
                                                    Dmin_local[j, m, 3] = i1
                                                    Dmin_local[j, m, 4] = i0
                                                    Dmin_local[j, m, 5] = i_1
                                                    Dmin_local[j, m, 6] = C

                                                xs_by_j_local[j].append(i_T)
                                                ys_by_j_local[j].append(E0)
                                                flag = 1

                                        if E0 > 0 and flag == 0:
                                            grey_segments_local.append([(i_T, E0), (i_T + 1, E0)])

        for j in range(1, 27):
            if xs_by_j[j]:
                plt.scatter(xs_by_j[j], ys_by_j[j], s=80, c=F[j], marker=".", linewidths=0)
        if grey_segments:
            lc = LineCollection(grey_segments, colors="#C0BCBC", linewidths=1)
            plt.gca().add_collection(lc)

        print("..............   plotting   .................")
        print("possible ET: ", i_T, "real ET: ", i_T1)
        print(
            "{0:10}{1:5}{2:20}{3:16}{4:8}{5:5}{6:5}{7:5}{8:5}{9:5}{10:5}{11:5}".format(
                "particle",
                "",
                "",
                "     theory: E",
                "   total ",
                "  i4",
                "  i3",
                "  i2",
                "  i1",
                "  i0",
                " i-1",
                "  C",
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
                # NOTE: left as-is (this is outside the hot scan loop)
                D_i_N[j] = float(i_Emax[j, 0]) * 100 / Cnt[m]
                Emax[j, m] = float(str(Emax[j, m])[:p_g])
                Emin[j, m] = float(str(Emin[j, m])[:p_g])
                E_mean = float(str(E_mean)[:p_g])

                print(
                    "{0:10}{1:5}{2:20}{3:16}{4:8}{5:5}{6:5}{7:5}{8:5}{9:5}{10:5}{11:5}".format(
                        "",
                        "max   ",
                        m_max,
                        Emax[j, m],
                        i_Emax[j, m],
                        Dmax[j, m, 0] / 2,
                        Dmax[j, m, 1] / 2,
                        Dmax[j, m, 2] / 2,
                        Dmax[j, m, 3] / 2,
                        Dmax[j, m, 4] / 2,
                        Dmax[j, m, 5] / 2,
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
                        "",
                        "min   ",
                        m_min,
                        Emin[j, m],
                        i_Emin[j, m],
                        Dmin[j, m, 0] / 2,
                        Dmin[j, m, 1] / 2,
                        Dmin[j, m, 2] / 2,
                        Dmin[j, m, 3] / 2,
                        Dmin[j, m, 4] / 2,
                        Dmin[j, m, 5] / 2,
                        Dmin[j, m, 6] / 2,
                    ),
                    file=f,
                )
                print(
                    "{0:10}{1:5}{2:20}{3:15}{4:8}{5:5}{6:8}{7:9}{8:10}{9:2}".format(
                        "",
                        "",
                        "",
                        "         ∆ abs(i)",
                        abs(Di_E),
                        "  Cts",
                        Cnt_,
                        "  ∆i/(2pi)",
                        D_i_c,
                        " %",
                    ),
                    file=f,
                )

                if flag == 1:
                    if Op == 1:
                        X = [0, 4, 7, 9, 6, 9, 5, 13, -21, -12, 5, -19, 6, -28, -17, -21, -13, -20, -12, -10, -37, -9, -5.5, -22, -11, 0, 0, 0]
                        Y = -20
                        fs = 12
                    if Op == 2:
                        X = [0, 1.2, 2.5, 3, 1, 1, 0.7, 1]
                        Y = -4
                        fs = 16
                    if Op == 3:
                        X = [0, 0.2, 0.2, 0.2]
                        Y = -1
                        fs = 16
                    if Op == 4:
                        X = [0, 2, 4, 5, 3, 4, 4, 8, -13, -7, 2, -11, 2, -14, -7, -11, -7, -12, -8, -3, -0, -0, -0, -22, -11, 0, 0, 0]
                        Y = -0
                        fs = 16
                    if Op == 5:
                        X = [0, 5, 10, 15, 20, 5, 9, 30, -35, -20, 7, -40, 25, 15, 4, 35, 20, -13, -10, -7, 20, 20, 10, -10, -10, -20, 0, 0, 0]
                        X[j] *= 2
                        Y = -50
                        fs = 12
                    plt.text(i_Emin[j, m] + 10000 * X[j], E_mean + Y, Obj[j][9], fontsize=fs, color=F[j])
                    flag = 2

            if j > 1 and m == 511 and i_Emax[j, 516] > 0:
                D_i_c = round(float(i_Emax[j, 0]) * 100 / i_Emax[j, 516], 5)
                D = str(D_i_c)
                D_i_c_[j] = "".rjust(10 - len(D), " ") + D + " %"
                print(
                    "{0:10}{1:5}{2:20}{3:15}{4:8}{5:5}{6:8}{7:9}{8:10}{9:2}".format(
                        "",
                        "total",
                        "",
                        "             Σ ∆i",
                        i_Emax[j, 0],
                        "  Cts",
                        i_Emax[j, 516],
                        "  ∆i/(2pi)",
                        D_i_c,
                        " %",
                        "\n",
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

            x_a = 0
            x_m = i_T * 1 / 5
            plt.ylabel("Energy in $m_e$")
            plt.xlabel("N")
            plt.xlim(-10000, i_T + 30000)

            i2_ = float(TWO_PI) ** 2
            i3_ = float(TWO_PI) ** 3
            i4_ = float(TWO_PI) ** 4
            i5_ = 1 / 2 * (i4_ + i3_ + i2_)
            i6_ = i4_ + i3_ + i2_
            i7_ = 5 / 2 * i4_ - 3 / 2 * i3_ - 1 / 2 * i3_
            i8_ = 3 / 2 * i4_ + i3_ + i3_
            i9_ = 2 * i4_ + 2 * i3_ + 2 * i3_
            i10_ = 3 / 2 * i4_ + 1 / 2 * i3_ - 1 / 2 * i3_

            if Op in [1]:
                plt.plot([x_a, x_m], [i4_, i4_], "k", linewidth=1)
                plt.text(x_a, i4_ + 15, r"$(2\pi)^4$", fontsize=12, color="blue")
            if Op in [5]:
                plt.plot([x_a, x_m], [i7_, i7_], "k", linewidth=1)
                plt.text(x_a, i7_ + 15, r"$5/2(2\pi)^4-3/2(2\pi)^3-1/2(2\pi)^2$", fontsize=12, color="blue")
            if Op in [5]:
                plt.plot([x_a, 2 * x_m], [i8_, i8_], "k", linewidth=1)
                plt.text(x_a, i8_ + 15, r"$3/2(2\pi)^4+(2\pi)^3+(2\pi)^2=", fontsize=12, color="blue")
            if Op in [5]:
                plt.plot([x_a, 3 * x_m], [i10_, i10_], "k", linewidth=1)
                plt.text(x_a, i10_ + 15, r"$3/2(2\pi)^4+1/2(2\pi)^3-1/2(2\pi)^2$", fontsize=12, color="blue")
            if Op in [1, 2]:
                plt.plot([x_a, x_m], [i3_, i3_], "k", linewidth=1)
                plt.text(x_a, i3_ + 15, r"$(2\pi)^3$", fontsize=12, color="blue")
            if Op in [1, 2, 3]:
                plt.plot([x_a, x_m], [i2_, i2_], "k", linewidth=1)
                plt.text(x_a, i2_ + 15, r"$(2\pi)^2$", fontsize=12, color="blue")
            if Op in [1]:
                plt.plot([x_a, x_m], [i5_, i5_], "k", linewidth=1)
                plt.text(x_a, i5_ + 15, r"$1/2((2\pi)^4+(2\pi)^3+(2\pi)^2)$", fontsize=12, color="blue")
            if Op in [1]:
                plt.plot([x_a, x_m], [i6_, i6_], "k", linewidth=1)
                plt.text(x_a, i6_ + 15, r"$(2\pi)^4+(2\pi)^3+(2\pi)^2$", fontsize=12, color="blue")

        # Legend blocks unchanged (same as before)
        if Op == 1:
            x_a = i_T * 0.65
            dx = i_T * 0.08
            i = -20
            for j in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 16, 17, 18, 19]:
                particle = str(Obj[j][9])
                plt.text(x_a, i, Obj[j][0])
                plt.text(x_a + dx, i, particle, color=F[j])
                plt.text(x_a + 2 * dx, i, i_Emax[j, 0])
                plt.text(x_a + 3 * dx, i, D_i_c_[j])
                i += 100
            plt.text(x_a + 2 * dx, i, "  ∆i  ")
            plt.text(x_a + 3 * dx, i, " ∆i/(2pi) ")

        if Op == 2:
            x_m = i_T * 1 / 5
            plt.xlim(-10000, i_T + 30000)
            x_a = i_T * 0.70
            dx = i_T * 0.10
            i = 0
            for j in [1, 2, 3, 4, 5, 6, 7]:
                particle = str(Obj[j][9])
                plt.text(x_a, i, Obj[j][0])
                plt.text(x_a + dx, i, particle, color=F[j])
                plt.text(x_a + 2 * dx, i, i_Emax[j, 0])
                plt.text(x_a + 3 * dx, i, D_i_c_[j])
                i += 20
            plt.text(x_a + 2 * dx, i, "  ∆i  ")
            plt.text(x_a + 3 * dx, i, " ∆i/(2pi) ")

        if Op == 3:
            x_m = i_T * 1 / 5
            plt.xlim(-1000, i_T + 10000)
            x_a = i_T * 0.89
            dx = i_T * 0.10
            i = 0
            for j in [1, 2, 3, 4, 5, 6, 7]:
                particle = str(Obj[j][9])
                plt.text(x_a, i, Obj[j][0])
                plt.text(x_a + dx, i, particle, color=F[j])

                D = D_i_c_[j]
                plt.text(x_a + 2 * dx, i, i_Emax[j, 0])
                plt.text(x_a + 3 * dx, i, D)

                print(D_i_c_[j])
                i += 5
            plt.text(x_a + 2 * dx, i, "  ∆i  ")
            plt.text(x_a + 3 * dx, i, " ∆i/(2pi) ")

        if Op == 4:
            x_m = i_T * 1 / 5
            plt.xlim(0, i_T)
            fs = 14
            x_a = i_T * 0.70
            dx = i_T * 0.10
            i = 0
            m_H = 1836.152673426 + 1
            plt.plot([1, i_T], [m_H, m_H], "k", linewidth=1)
            plt.text(1, 1837, "$m_{Proton} + m_e$", fontsize=fs, color="blue")
            for j in [17, 19]:
                plt.text(1050, float(Obj[j][2]), Obj[j][0], fontsize=fs, color=F[j])
                particle = str(Obj[j][9])
                plt.text(x_a, i, Obj[j][0])
                plt.text(x_a + dx, i, particle, color=F[j])
                plt.text(x_a + 2 * dx, i, i_Emax[j, 0])
                plt.text(x_a + 3 * dx, i, D_i_c_[j])
                i += 20

        if Op == 5:
            x_a = i_T * 0.65
            dx = i_T * 0.085
            i = 1500
            for j in [15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25]:
                particle = str(Obj[j][9])
                if j == 26:
                    particle = ""
                plt.text(x_a, i, Obj[j][0])
                plt.text(x_a + dx, i, particle, color=F[j])
                plt.text(x_a + 2 * dx, i, i_Emax[j, 0])
                plt.text(x_a + 3 * dx, i, D_i_c_[j])
                i += 100
            plt.text(x_a + 2 * dx, i, "  ∆i  ")
            plt.text(x_a + 3 * dx, i, " ∆i/(2pi) ")

    fig = plt.gcf()
    fig.set_size_inches(10, 6)
    
    fig.savefig(os.path.join(RESULTS_DIR, f"{config.name}.png"), dpi=100)
    
    end_time_stamp = datetime.datetime.now()
    delta = end_time_stamp - start_time_stamp
    print(f"took {delta.seconds} seconds")
    
    if not no_show:
        plt.show()


if __name__ == "__main__":
    main()