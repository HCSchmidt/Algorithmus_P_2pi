import cmath
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

from ..cli import ScanSector
from ..constants import get_colors
  

def draw_points(xs_by_j, ys_by_j, grey_segments, N_segments, Obj, Charge, Ch, show_H, show_N):      # Charge[0][Ch]                   #    um Charge ergänzt

    colors = get_colors()
    for j in range(1, 29):
        if Obj[j][6] == Charge[0][Ch] or " " == Charge[0][Ch] and not j == show_H and not j == show_N:  
            print(Obj[j][6])
            if xs_by_j[j]:    
               plt.scatter(xs_by_j[j], ys_by_j[j], s=40, c=colors[j], marker=".", linewidths=0)
               
    if grey_segments:
        lc = LineCollection(grey_segments, colors="#6F6E6E", linewidths=2)
        plt.gca().add_collection(lc)

#    if N_segments:
#        lc = LineCollection(N_segments, colors="#0F0F0F", linewidths=2)
#        plt.gca().add_collection(lc)
 
 
    
                
def add_reference_lines(sector: ScanSector, i_T: int):
    x_a = 0
    x_m = i_T * 1 / 5

    plt.ylabel("Energy in $m_e$")
    plt.xlabel("N")

    TWO_PI = float(2 * cmath.pi)
    i2_ = TWO_PI ** 2
    i1_ = TWO_PI 
    i3_ = TWO_PI ** 3
    i4_ = TWO_PI ** 4
    i5_ = 0.5 * (i4_ + i3_ + i2_)
    i6_ = i4_ + i3_ + i2_
    i7_= 2.5* i4_ - 1.5 * i3_ - 0.5 * i2_;
    i9_ = 2 * i4_ + 2 * i3_ + 1.5 * i2_;
    i10_ = 1.5 * i4_ + 0.5 * i3_ + 0.5 * i2_;
    i11_ = 2 * i4_ - 2 * i3_  - 2 * i2_ - 2 * i1_ 
    i12_ = 1 * i4_ - 2 * i3_  - 2 * i2_ - 2 * i1_ 
    i13_ = 3/2 * i4_ - 2 * i3_  - 2 * i2_ - 2 * i1_   
    i14_ = 1/2 * i4_ - 2 * i3_  - 2 * i2_ - 2 * i1_ 

    if sector is ScanSector.broad:
        plt.plot([x_a, x_m], [i4_, i4_], "k", linewidth=1)
        plt.text(x_a, i4_ + 15, r"$(2\pi)^4$", fontsize=12, color="blue")
        plt.plot([x_a, x_m], [i3_, i3_], "k", linewidth=1)
        plt.text(x_a, i3_ + 15, r"$(2\pi)^3$", fontsize=12, color="blue")
        plt.plot([x_a, x_m], [i5_, i5_], "k", linewidth=1)
        plt.text(x_a, i5_ + 15, r"$1/2((2\pi)^4+(2\pi)^3+(2\pi)^2)$", fontsize=12, color="blue")
        plt.plot([x_a, x_m], [i6_, i6_], "k", linewidth=1)
        plt.text(x_a, i6_ + 15, r"$(2\pi)^4+(2\pi)^3+(2\pi)^2$", fontsize=12, color="blue")

    if sector is ScanSector.E112P:
        plt.plot([x_a, x_m], [i6_, i6_], "k", linewidth=1)
        plt.text(x_a, i6_ + 15, r"$(2\pi)^4+(2\pi)^3+(2\pi)^2$", fontsize=12, color="blue")

    if sector is ScanSector.heavy:
        plt.plot([x_a, 4.5 * x_m], [i7_, i7_], "k",linewidth=1); 
        plt.text(x_a, i7_ + 15, r"$5/2(2\pi)^4-2(2\pi)^3-2(2\pi)^2-2(2\pi)$", fontsize=12, color='blue')
        plt.plot([x_a, 3.8 * x_m], [i11_, i11_], "k", linewidth=1); 
        plt.text(x_a, i11_ + 15, r"$2(2\pi)^4-2(2\pi)^3-2(2\pi)^2-2(2\pi)$", fontsize=12, color='blue')
        plt.plot([1.8 * x_m, 5 * x_m], [i12_, i12_], "k", linewidth=1); 
        plt.text(3 * x_m, i12_ + 15, r"$(2\pi)^4-2(2\pi)^3-2(2\pi)^2-2(2\pi)$", fontsize=12, color='blue')
        plt.plot([2.8 * x_m, 5 * x_m], [i13_, i13_], "k", linewidth=1); 
        plt.text(3.3 * x_m, i13_ - 150, r"$3/2(2\pi)^4-2(2\pi)^3-2(2\pi)^2-2(2\pi)$", fontsize=12, color='blue')
        plt.plot([0.8 * x_m, 5 * x_m], [i14_, i14_], "k", linewidth=1); 
        plt.text(3 * x_m, i14_ + 15, r"$1/2(2\pi)^4-2(2\pi)^3-2(2\pi)^2-2(2\pi)$", fontsize=12, color='blue')

    if sector is ScanSector.E222:
        plt.plot([x_a, 4.5 * x_m], [i11_, i11_], "k", linewidth=1); 
        plt.text(x_a, i11_ + 15, r"$2(2\pi)^4-2(2\pi)^3-2(2\pi)^2-2(2\pi)$", fontsize=12, color='blue')
        plt.plot([x_a, 2.0 * x_m], [i13_, i13_], "k", linewidth=1); 
        plt.text(x_a, i13_ + 15, r"$3/2(2\pi)^4-2(2\pi)^3-2(2\pi)^2-2(2\pi)$", fontsize=12, color='blue')
        plt.plot([2.5 * x_m, 5 * x_m], [i12_, i12_], "k", linewidth=1); 
        plt.text(3 * x_m, i12_ + 15, r"$(2\pi)^4-2(2\pi)^3-2(2\pi)^2-2(2\pi)$", fontsize=12, color='blue')
        plt.plot([1.5 * x_m, 5 * x_m], [i14_, i14_], "k", linewidth=1); 
        plt.text(3 * x_m, i14_ + 15, r"$1/2(2\pi)^4-2(2\pi)^3-2(2\pi)^2-2(2\pi)$", fontsize=12, color='blue')

    if sector is ScanSector.E333D:
        plt.plot([x_a, x_m], [i6_, i6_], "k", linewidth=1)
        plt.text(x_a, i6_ + 15, r"$(2\pi)^4+(2\pi)^3+(2\pi)^2$", fontsize=12, color="blue")      
    
    if sector in (ScanSector.light):
        plt.plot([x_a, x_m], [i3_, i3_], "k", linewidth=1)
        plt.text(x_a, i3_ + 15, r"$(2\pi)^3$", fontsize=12, color="blue")
        plt.plot([x_a, x_m], [i2_, i2_], "k", linewidth=1)
        plt.text(x_a, i2_ + 15, r"$(2\pi)^2$", fontsize=12, color="blue")

    if sector in (ScanSector.minimal):
        plt.plot([x_a, x_m], [i1_, i1_], "k", linewidth=1)
        plt.text(x_a, i1_ +1 , r"$(2\pi)$", fontsize=12, color="blue")

    # x-limits similar to your previous plots
    if sector is ScanSector.minimal:
        plt.xlim(-1000, i_T + 10000)
    elif sector is ScanSector.nucleon:
        plt.xlim(0, i_T)
    else:
        plt.xlim(-10000, i_T + 30000)

def label_offsets_for_sector(sector: ScanSector, j: int):
    """Return (X, Y, fs) where X is a whole list so we can use X[j]."""
    if sector is ScanSector.broad:
        X = [0, 4, 7, 9,8 , 9, 5, 13, -21, -12, 7, -19, 8, -28, -17, -21, -13, -20, -12, -10, -37, -9, -5.5, -22, -11, 0, 0, 0]
        return X, -20, 12

    if sector is ScanSector.E112P:                                          #H   
        X = [0, 4, 7, 9, 8, 9, 5, 13, -21, -12, 7, -19, 8, -28, -17, -4, 6, 1, -4, -4, -5, 0, 0, 0, 0, 0, 0, 0]
        return X, 0, 16
  
#    if sector is ScanSector.E122:       # test  schräg                                  #H   
#        X = [0, 4, 7, 9, 8, 9, 5, 13, -21, -12, 7, -19, 8, -28, -17, -4, 6, 1, -4, -4, -5, 0, 0, 0, 0, 0, 0, 0]
#        return X, 0, 16

    if sector is ScanSector.light:
        X = [0, 1.2, 2.5, 3,-3, 1, -4, 1]
        return X, -4, 16

    if sector is ScanSector.minimal:
        X = [0, 0.1, 0.3, 0.1]
        return X, -1, 16

    if sector is ScanSector.nucleon:
        X = [0, 2, 4, 5, 3, 4, 4, 8, -13, -7, 2, -11, 2, -14, -7, 1, 1, 1, 1, 1, 1, -0, -0, -22, -11, 0, 0, 0]
        return X, -0, 16

    if sector is ScanSector.heavy:                                                  #H
        X = [0, 5, 10, 15, 20, 24, 20, 38, -35, -20, 18, -38, -15, 5, 22, 25, -20, -30, -19, -14, 20, 20, 10, -10, -30, -15, -5, 0, 0]
        X[j] *= 2
        return X, -30, 12
    
    if sector is ScanSector.E222:                                                  #H
        X = [0, 5, 10, 15, 20, 24, 20, 38, -35, -20, 18, -38, -15, 5, 22, 25, -20, -30, -19, -14, 20, 20, 10, -10, -30, -15, -5, 0, 0]
        X[j] *= 2
        return X, -30, 12
    
    if sector is ScanSector.E333U:                                          #P                       #tau
        X = [0, 5, 10, 15, 20, 5, 9, 30, -35, -20, 7, -40, 25, 0, 0, 60, 40, 3 ,12, 25, 35, 50, 10, -10, -40, -25, -9, -15, -15, 0, 0]
        X[j] *= 4
        return X, -50, 12
    
    if sector is ScanSector.E333D:                                          #P                       #tau
        X = [0, 4, 6, 8, -3, 4, 4, 20, -10, 4, 23, -19, -10, 4, 16, -18, -9, -24, -18, -15, -8, -15, -0, -22, -11, 0, 0, 0]
        X[j] *= 5
        return X, -0, 16  

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

"""

u: 4.18(-0.51)(0.96) &  $\overline{E} \approx 2/3*2 = 4.2 $
1/3*2+ 1/3 + 2*2^(-1) = 2 $ \\
2/3*2+ 2/3 + 2^(-1)= 2.5 $ \ \ besser 3.67  \\
2+ 1 + 2^(-1) = 3.5 \ = \ (m_n - m_p) / (m_p) \cdot 608.8 s / 86164 s $ \ \ Zerfall von n \\
\hline
d: 9.14(-0.33)(0.94)  & $\overline{E} \approx 3/2*2 = 9.4$ 
 4/3 + 1/3*2^(-1) = 1.5 $ \\
3/2*2+2/3+ 1/3*2^(-1)=3.833$\\
3/2*2+2/3+ 0/3*2^(-1) = 3.66$\\
2*2+2+3*2^(-1) = 7.5 $\\
1 + 3/2*2^(-1) = 1.75$\\
\hline
s: 182.8(-6.6)(16.8) & $\overline{E} \approx 2/3*2^3 + 1/2*2^2 = 185.1 $ \ \ \ $\Delta E \approx 1/3*2^2 + 4/3 *2^1+4/3+4/3*2^(-1)+... = 23.1$\\

176.152/3*2^3+1/3*2^2+1/3*2^1+1/3+1/3*2^(-1)=7.833 \\
2/3*2^3+2/3*2^2+*2^1+1+*2^(-1)=17.5 $\\
2^3+2^2+2^1+3+2^(-1)= 17.5$\\
1/2*2^2 + 1/2 *2^1+1/2+1/2*2^(-1) =23.4$\\
\hline
c: 2485(-39)(39) & $\overline{E} \approx 3/2*2^4+1/2*2^3 = 188.6$
\ \ \ $\Delta E \approx 2*2^2 + 1 = 78$ \\
3/2*2^4+1/2*2^3+1/2*2^2+1/2*2^1+1/2+1/2*2^(-1)= 31.75 $\\
3/2*2^4+1/2*2^3+3/2*2^2+1/2*2^1+1/2+1/2*2^(-1)=35.75 $\\
2*2^4+2*2^3+2*2^2+2*2^1+5+2*2^(-1) = 66 $\\
 2*2^2 - 1 = 78$ \\
\hline
b: 8100(-60)(80)  & $\overline{E} \approx 3/4*2^5+1/2*2^4 = 8124$ \ \ \ $\Delta E \approx 2/3*2^3 + 1/2*2^2 + *2^1+... = 140$ \\
&  
3/4*2^5+1/2*2^4+1/3*2^3+0/2*2^2+0/3*2^1+1+...= 35.66 $\\
3/4*2^5+1/2*2^4+1/3*2^3+1/2*2^2+1*2^1+0+...= 38.66 $\\
2^5+2^4+0*2^3+*2^2+2*2^1+2 = 58$ \\
2^3 + 2*2^2+ 4*2^1+3= 140$ \\
\hline
t: 337710(570) & $\overline{E} \approx 3/4*2^7+2/3*2^6+2/3*2^5= 337496$ \ \ \ $\Delta E \approx 2^4 + 2*2^3+2*2^2+... = 2 \cdot 570$\\
3/4*2^7+2/3*2^6+2/3*2^5+1/3*2^4+2/3*2^3+0*2^2+0*2= 170.66 $\\
3/4*2^7+2/3*2^6+2/3*2^5+2/3*2^4+2^3+1/3*2^2+2= 182$\\
2^7+2^6+2^5+2*2^4+2^3+2^2+2*2^1+1 = 273 $\\
2^4 + 2*2^3+2*2^2 = 1140$\\
\hline
Higgs: 244830(210)  & $\overline{E} \approx 4*2^6+0*2^5+3/4*2^4=244947$ \ \ \ \ $\Delta E \approx 2*2^3 + 2*2^2 + 1/2*2^1 = 2 \cdot 210$\\
& 
4*2^6+0*2^5+3/4*2^4+4/3*2^3+0*2^2+3/4*2^1+1 = 281.16 $\\
& 
4*2^6 + 0*2^5+3/4*2^4+1/4*2^3+3/4*2^2+0*2^1+1 =274 $\\
\hline

K$^{++}$: 966.102(21) & $E_{K} = *2^4- E_n=  968.60 $ \\
K$^0$: 973.800(26) \\ 

K*$^{+-}$: 1745.2(1) & $E_{K^*} = 3/2*2^4- E_n = 1747.87 $ \\
K*$^0$: 1752.6(1)\\ 

"""
