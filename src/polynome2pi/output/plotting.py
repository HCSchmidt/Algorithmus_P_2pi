import cmath
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

from ..cli import ScanSector
from ..constants import get_colors
from ..particles import Particle


# ---------------------------------------------------------------------
# Scatter + grey background segments
# ---------------------------------------------------------------------

def draw_points(xs_by_particle, ys_by_particle, grey_segments):
    """
    xs_by_particle / ys_by_particle:
        dict[str, list[float]] keyed by particle.key
    """
    colors = get_colors()

    for idx, (key, xs) in enumerate(xs_by_particle.items(), start=1):
        if xs:
            plt.scatter(
                xs,
                ys_by_particle[key],
                s=40,
                c=colors[idx],
                marker=".",
                linewidths=0,
            )

    if grey_segments:
        lc = LineCollection(grey_segments, colors="#6F6E6E", linewidths=2)
        plt.gca().add_collection(lc)


# ---------------------------------------------------------------------
# Reference energy lines
# ---------------------------------------------------------------------

def add_reference_lines(sector: ScanSector, i_T: int):
    x_a = 0
    x_m = i_T / 5

    plt.ylabel("Energy in $m_e$")
    plt.xlabel("N")

    TWO_PI = float(2 * cmath.pi)
    i2_ = TWO_PI**2
    i3_ = TWO_PI**3
    i4_ = TWO_PI**4
    i5_ = 0.5 * (i4_ + i3_ + i2_)
    i6_ = i4_ + i3_ + i2_
    i7_ = 2.5 * i4_ - 1.5 * i3_ - 0.5 * i2_
    i9_ = 2 * i4_ + 2 * i3_ + 1.5 * i2_
    i10_ = 1.5 * i4_ + 0.5 * i3_ + 0.5 * i2_

    if sector is ScanSector.broad:
        _line(x_a, x_m, i4_, "$(2\\pi)^4$")
        _line(x_a, x_m, i3_, "$(2\\pi)^3$")
        _line(x_a, x_m, i2_, "$(2\\pi)^2$")
        _line(x_a, x_m, i5_, "$1/2((2\\pi)^4+(2\\pi)^3+(2\\pi)^2)$")
        _line(x_a, x_m, i6_, "$(2\\pi)^4+(2\\pi)^3+(2\\pi)^2$")

    if sector is ScanSector.heavy:
        _line(x_a, 5 * x_m, i7_, "$5/2(2\\pi)^4-3/2(2\\pi)^3-1/2(2\\pi)^2$")
        _line(x_a, 4 * x_m, i9_, "$2(2\\pi)^4+2(2\\pi)^3+3/2(2\\pi)^2$")
        _line(x_a, 3 * x_m, i10_, "$3/2(2\\pi)^4+1/2(2\\pi)^3+1/2(2\\pi)^2$")

    if sector in (ScanSector.broad, ScanSector.light):
        _line(x_a, x_m, i3_, "$(2\\pi)^3$")

    if sector in (ScanSector.broad, ScanSector.light, ScanSector.minimal):
        _line(x_a, x_m, i2_, "$(2\\pi)^2$")

    if sector is ScanSector.minimal:
        plt.xlim(-1000, i_T + 10000)
    elif sector is ScanSector.nucleon:
        plt.xlim(0, i_T)
    else:
        plt.xlim(-10000, i_T + 30000)


def _line(x0, x1, y, label):
    plt.plot([x0, x1], [y, y], "k", linewidth=1)
    plt.text(x0, y + 15, label, fontsize=12, color="blue")


# ---------------------------------------------------------------------
# Legend panels
# ---------------------------------------------------------------------

def add_legend_panels(
    sector: ScanSector,
    i_T: int,
    particles: dict[str, Particle],
    i_Emax: dict[str, int],
    D_i_c_: dict[str, str],
):
    colors = get_colors()
    items = list(particles.values())

    if sector is ScanSector.broad:
        _legend_block(i_T, items, colors, i_Emax, D_i_c_, y_start=-20, dy=100)

    if sector is ScanSector.light:
        _legend_block(i_T, items[:7], colors, i_Emax, D_i_c_, y_start=0, dy=20)

    if sector is ScanSector.minimal:
        _legend_block(i_T, items[:7], colors, i_Emax, D_i_c_, y_start=0, dy=5)

    if sector is ScanSector.nucleon:
        _legend_nucleon(i_T, particles, colors, i_Emax, D_i_c_)

    if sector is ScanSector.heavy:
        _legend_block(i_T, items[14:25], colors, i_Emax, D_i_c_, y_start=1500, dy=100)


def _legend_block(i_T, items, colors, i_Emax, D_i_c_, y_start, dy):
    x_a = i_T * 0.65
    dx = i_T * 0.085
    y = y_start

    for idx, p in enumerate(items, start=1):
        plt.text(x_a, y, p.name)
        plt.text(x_a + dx, y, p.symbol, color=colors[idx])
        plt.text(x_a + 2 * dx, y, i_Emax.get(p.key, 0))
        plt.text(x_a + 3 * dx, y, D_i_c_.get(p.key, ""))
        y += dy

    plt.text(x_a + 2 * dx, y, "∆i")
    plt.text(x_a + 3 * dx, y, "∆i/(2π)")


def _legend_nucleon(i_T, particles, colors, i_Emax, D_i_c_):
    fs = 14
    x_a = i_T * 0.7
    dx = i_T * 0.1
    y = 0

    proton = particles["proton"]
    neutron = particles["neutron"]

    m_H = proton.theory_E + 1
    plt.plot([1, i_T], [m_H, m_H], "k", linewidth=1)
    plt.text(1, m_H + 1, "$m_{Proton} + m_e$", fontsize=fs, color="blue")

    for idx, p in enumerate((proton, neutron), start=1):
        plt.text(x_a, y, p.name)
        plt.text(x_a + dx, y, p.symbol, color=colors[idx])
        plt.text(x_a + 2 * dx, y, i_Emax.get(p.key, 0))
        plt.text(x_a + 3 * dx, y, D_i_c_.get(p.key, ""))
        y += 20


# ---------------------------------------------------------------------
# Particle labels in plot area
# ---------------------------------------------------------------------

def add_particle_labels(labels):
    """
    labels: iterable of tuples:
        (sector, particle, x_base, y)
    """
    if not labels:
        return

    colors = get_colors()

    # stable mapping from particle key -> color index
    # (uses insertion order of PARTICLES dict)
    from ..particles import PARTICLES
    key_to_idx = {k: i for i, k in enumerate(PARTICLES.keys(), start=1)}

    for sector, particle, x_base, y in labels:
        X, Y, fs = label_offsets_for_sector(sector, particle.key)
        idx = key_to_idx.get(particle.key, 1)
        color = colors[idx]

        plt.text(
            x_base + 10000 * X,
            y + Y,
            particle.symbol,
            fontsize=fs,
            color=color,
        )


def label_offsets_for_sector(sector: ScanSector, particle_key: str):
    """
    Return (X_offset, Y_offset, font_size)
    """
    if sector is ScanSector.broad:
        return -20, -20, 12
    if sector is ScanSector.light:
        return 2, -4, 16
    if sector is ScanSector.minimal:
        return 0.2, -1, 16
    if sector is ScanSector.nucleon:
        return 2, 0, 16
    if sector is ScanSector.heavy:
        return 10, -50, 12

    raise ValueError(f"Unhandled sector: {sector}")