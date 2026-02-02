import cmath

import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

from .cli import ScanSector
from .constants import get_colors
    

def draw_points(xs_by_j, ys_by_j, grey_segments):
    colors = get_colors()
    for j in range(1, 27):
        if xs_by_j[j]:
            plt.scatter(xs_by_j[j], ys_by_j[j], s=80, c=colors[j], marker=".", linewidths=0)

    if grey_segments:
        lc = LineCollection(grey_segments, colors="#C0BCBC", linewidths=1)
        plt.gca().add_collection(lc)




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


def add_legend_panels(sector: ScanSector, i_T: int, Obj, i_Emax, D_i_c_):
    colors = get_colors()
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