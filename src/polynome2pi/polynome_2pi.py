import os
import datetime
from pathlib import Path
import cmath

import numpy as np

import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

from .cli import ScanSector, PolynomeConfig, select_preset_by_sector, parse_args
from .engine import run_scan
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

def main(argv=None):
    start_time_stamp = datetime.datetime.now()

    args = parse_args(argv)
    sector = args.sector
    config = select_preset_by_sector(sector)
    no_show = args.no_show

    # data
    Obj, obj_E, obj_min, obj_max = build_particle_table()
    D_i_N, D_i_c_, Cnt, Emax, Emin, i_Emax, i_Emin, Dmax, Dmin = allocate_result_arrays()
    xs_by_j, ys_by_j, grey_segments = init_plot_buffers()

    # scan
    i_T, i_T1, _mmax = run_scan(
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