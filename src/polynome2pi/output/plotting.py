import cmath

import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

from ..cli import ScanSector
from ..constants import get_colors
    
def draw_points(xs_by_j, ys_by_j, grey_segments,Obj, Charge, Ch, show_H):      # Charge[0][Ch]                   #    um Charge ergänzt
    colors = get_colors()
    for j in range(1, 27):
        if Obj[j][6] == Charge[0][Ch] or " " == Charge[0][Ch] and (not j == show_H[0] ): #   and  not j == 18):  
            print(Obj[j][6])
            if xs_by_j[j]:    
                plt.scatter(xs_by_j[j], ys_by_j[j], s=40, c=colors[j], marker=".", linewidths=0)


    if grey_segments:
        lc = LineCollection(grey_segments, colors="#6F6E6E", linewidths=2)
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
    i7_= 2.5* i4_ - 1.5 * i3_ - 0.5 * i2_;
    i9_ = 2 * i4_ + 2 * i3_ + 1.5 * i2_;
    i10_ = 1.5 * i4_ + 0.5 * i3_ + 0.5 * i2_;

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


    if sector is ScanSector.E112P:
    #    plt.plot([x_a, x_m], [i4_, i4_], "k", linewidth=1)
    #    plt.text(x_a, i4_ + 15, r"$(2\pi)^4$", fontsize=12, color="blue")
        plt.plot([x_a, x_m], [i6_, i6_], "k", linewidth=1)
        plt.text(x_a, i6_ + 15, r"$(2\pi)^4+(2\pi)^3+(2\pi)^2$", fontsize=12, color="blue")


    if sector is ScanSector.heavy:
        plt.plot([x_a, 5 * x_m], [i7_, i7_], "k", linewidth=1); 
        plt.text(x_a, i7_+ 15, r"$5/2(2\pi)^4-3/2(2\pi)^3-1/2(2\pi)^2$", fontsize=12, color='blue')
        plt.plot([x_a, 3.2 * x_m], [i9_, i9_], "k",linewidth=1); 
        plt.text(x_a, i9_ + 15, r"$2(2\pi)^4+2(2\pi)^3+3/2(2\pi)^2$", fontsize=12, color='blue')
        plt.plot([x_a, 3 * x_m], [i10_, i10_], "k", linewidth=1); 
        plt.text(x_a, i10_ + 15, r"$3/2(2\pi)^4+1/2(2\pi)^3+1/2(2\pi)^2$", fontsize=12, color='blue')

    
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
        i = -80
        for jj in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 16, 17, 18, 19]:
            particle = str(Obj[jj][9])
            plt.text(x_a, i, Obj[jj][0])
            plt.text(x_a + dx, i, particle, color=colors[jj])
            plt.text(x_a + 2 * dx, i, i_Emax[jj, 0])
            plt.text(x_a + 3 * dx, i, D_i_c_[jj])
            i += 100
        plt.text(x_a + 2 * dx, i, "  ∆i  ")
        plt.text(x_a + 3 * dx, i, " ∆i/(2pi) ")


    if sector is ScanSector.E112P:
        x_a = i_T * 0.5
        dx = i_T * 0.15
        pi = cmath.pi ; E1700= (2 * pi)**4 +(2 * pi)**2+ (2 * pi)**1
        i = E1700
        for jj in [ 15, 16, 17, 18, 19, 20]:
            particle = str(Obj[jj][9])
            plt.text(x_a, i, Obj[jj][0])
            plt.text(x_a + dx, i, particle, color=colors[jj])
            plt.text(x_a + 2 * dx, i, i_Emax[jj, 0])
            plt.text(x_a + 3 * dx, i, D_i_c_[jj])
            i += 15
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
        x_a = i_T * 0.62
        dx = i_T * 0.085
        i = 1410
        for jj in [15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26]:
            particle = str(Obj[jj][9])
            plt.text(x_a, i, Obj[jj][0])
            plt.text(x_a + dx, i, particle, color=colors[jj])
            plt.text(x_a + 2 * dx, i, i_Emax[jj, 0])
            plt.text(x_a + 3 * dx, i, D_i_c_[jj])
            i += 100
        plt.text(x_a + 2 * dx, i, "  ∆i  ")
        plt.text(x_a + 3 * dx, i, " ∆i/(2pi) ")
        
        
def label_offsets_for_sector(sector: ScanSector, j: int):
    """Return (X, Y, fs) where X is a whole list so we can use X[j]."""
    if sector is ScanSector.broad:
        X = [0, 4, 7, 9, 6, 9, 5, 13, -21, -12, 5, -19, 6, -28, -17, -21, -13, -20, -12, -10, -37, -9, -5.5, -22, -11, 0, 0, 0]
        return X, -20, 12

    if sector is ScanSector.E112P:
        X = [0, 4, 7, 9, 6, 9, 5, 13, -21, -12, 5, -19, 6, -28, -17, -6, 4, 1, -4, -4, -5, 0, 0, 0, 0, 0, 0, 0]
        return X, 0, 16

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
        X = [0, 5, 10, 15, 20, 5, 9, 30, -35, -20, 7, -40, 25, 15, 4, 35, 20, -13, -10, -7, 20, 20, 10, -10, -15, -25, -5, 0, 0]
        X[j] *= 2
        return X, -50, 12

    raise ValueError(f"Unhandled sector: {sector}")

def add_particle_labels(labels):
    """Draw all particle labels returned by report.write_report().

    Expected `labels` format: iterable of dicts with keys:
      - sector, j, x_base, y, text, color
    """
    if not labels:
        return

    for item in labels:
        sector, j, x_base, y, text, color = item[:6]
        X, Y, fs = label_offsets_for_sector(sector, j)
        plt.text(
            x_base + 10000 * X[j],
            y + Y,
            text,
            fontsize=fs,
            color=color,
        )
        
