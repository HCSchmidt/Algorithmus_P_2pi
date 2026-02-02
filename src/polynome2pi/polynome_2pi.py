import os
import datetime
from pathlib import Path
import cmath

import numpy as np

import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

from .cli import ScanSector, PolynomeConfig, select_preset_by_sector, parse_args
from .engine import PolynomeEngine
from .constants import build_particle_table
from .plotting import add_reference_lines, draw_points, add_legend_panels
from .report import write_report_and_labels
RESULTS_DIR = Path("results")
RESULTS_DIR.mkdir(exist_ok=True)



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
    sector: ScanSector,
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
    
    is_heavy = (sector is ScanSector.heavy)
    is_nucleon = (sector is ScanSector.nucleon)

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









# ------------------------- main orchestration -------------------------

def main(argv=None):
    start_time_stamp = datetime.datetime.now()

    args = parse_args(argv)
    sector = args.sector
    config = select_preset_by_sector(sector)
    no_show = args.no_show

    engine = PolynomeEngine()

    # data
    Obj, obj_E, obj_min, obj_max = build_particle_table()
    D_i_N, D_i_c_, Cnt, Emax, Emin, i_Emax, i_Emin, Dmax, Dmin = allocate_result_arrays()
    xs_by_j, ys_by_j, grey_segments = init_plot_buffers()

    # scan
    i_T, i_T1, _mmax = run_scan(
        engine=engine,
        config=config,
        sector=sector,
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
    draw_points(xs_by_j, ys_by_j, grey_segments)

    # report + particle labels
    print("possible ET:", i_T, "real ET:", i_T1)
    write_report_and_labels(
        file_path = RESULTS_DIR / f"{config.name}.txt",
        sector=sector,
        Obj=Obj,
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
    add_legend_panels(sector, i_T, Obj, i_Emax, D_i_c_)

    # save/show
    fig = plt.gcf()
    fig.set_size_inches(10, 6)
    fig.savefig(RESULTS_DIR / f"{config.name}.png", dpi=100)

    end_time_stamp = datetime.datetime.now()
    delta = end_time_stamp - start_time_stamp
    print(f"took {delta.seconds} seconds")

    if not no_show:
        plt.show()


if __name__ == "__main__":
    main()