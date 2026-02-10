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
from .output.plotting import add_reference_lines, draw_points, add_legend_panels, add_particle_labels
from .output.report import write_report
from .output.config import RESULTS_DIR
from .initialize import init_result_arrays, init_plot_buffers
from .utils import open_image

Path(RESULTS_DIR).mkdir(exist_ok=True)

def main(argv=None):
    start_time_stamp = datetime.datetime.now()

    args = parse_args(argv)
    sector = args.sector
    config = select_preset_by_sector(sector)
    no_show = args.no_show
    
    Ch = 6                                                  #    für die Ladung evtl. zu gebrauchen
    # parse_args(argvc).Charge
    # = argvc; # argsc.Charge         
    # data

    Obj, obj_E, obj_min, obj_max = build_particle_table()
    D_i_N, D_i_c_, Cnt, Emax, Emin, i_Emax, i_Emin, Dmax, Dmin = init_result_arrays()
    xs_by_j, ys_by_j, grey_segments = init_plot_buffers()

    Charge = [
        ["1","2/3","1/3","0","-1/3","-2/3","-1","+-1"],       #für die Ladung 
        [" 1 ","2 3","1 3"," 0 ","-1 3","-2 3","-1 ","+-1"],  # den Dateiname vom plot 
       ]  

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
    draw_points(xs_by_j, ys_by_j, grey_segments,Obj, Charge, Ch)      #  ergänzt um Charge

    # report + particle labels
    print("possible ET:", i_T, "real ET:", i_T1)
    labels = write_report(
        config.name,
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
    add_particle_labels(labels)

    # save/show
    # output_png = os.path.join(RESULTS_DIR, f"{config.name}.png")
    output_png = os.path.join(RESULTS_DIR, f"{config.name+ " "+ Charge[1][Ch] + " "}.png")   # ergänzt um Charge

    fig = plt.gcf()
    fig.set_size_inches(10, 6)
    fig.savefig(output_png, dpi=100)

    end_time_stamp = datetime.datetime.now()
    delta = end_time_stamp - start_time_stamp
    print(f"took {delta.seconds} seconds")

    if not no_show:
        open_image(output_png)

if __name__ == "__main__":

    main()

