import os
import datetime
from pathlib import Path
import matplotlib.pyplot as plt

from .cli import ScanSector, PolynomeConfig, select_preset_by_sector, parse_args
from .engine import run_scan
from .constants import build_particle_table
from .output.plotting import add_reference_lines, draw_points, add_particle_labels
from .output.report import write_report
from .output.config import RESULTS_DIR
from .initialize import init_result_arrays, init_plot_buffers
from .utils import open_image
from .output.legend import add_legend_panels

Path(RESULTS_DIR).mkdir(exist_ok=True)

def main(argv=None, argC=None, argH=None, argV=None):
    start_time_stamp = datetime.datetime.now()

    args = parse_args(argv)
    sector = args.sector
    config = select_preset_by_sector(sector)
    no_show = args.no_show

    Charge = [
        [" ","1","2/3","1/3","0","-1/3","-2/3","-1","+-1"],       #für die Ladung 
        [" "," 1 ","2 3","1 3"," 0 ","-1 3","-2 3","-1 ","+-1"," "]    # den Dateiname vom plot 
    ]                                       
    for c in range(0,8):
        if Charge[0][c] == argC[1]: Ch = c

    show_H = -1; show_N = -1;                # H = argH[1]   # Atom with or without electrons
    if not argH[1] == "H": show_H = 18
    if not argH[2] == "N": show_N = 19       # and not j == show_H and not j == show_N
  
    Obj, obj_E, obj_min, obj_max = build_particle_table()
    D_i_N, D_i_c_, Cnt, Emax, Emin, i_Emax, i_Emin, Dmax, Dmin, PEmax, PEmin = init_result_arrays()
    xs_by_j, ys_by_j, grey_segments, N_segments, data = init_plot_buffers() 

    i_T, i_T1, mmax = run_scan(
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
        PEmax=PEmax,
        PEmin=PEmin,
        xs_by_j=xs_by_j,
        ys_by_j=ys_by_j,
        grey_segments=grey_segments,
        N_segments=N_segments,
        data=data
    )

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
        PEmax=PEmax,
        PEmin=PEmin,
    )

    # plot points (batched)
    draw_points(xs_by_j, ys_by_j, grey_segments, N_segments, Obj, Charge, Ch, show_H, show_N)      #  ergänzt um Charge
    add_reference_lines(sector, i_T)
    add_particle_labels(labels)

    output_png = os.path.join(RESULTS_DIR, f"{config.name+ "_"+ Charge[1][Ch] + "_"+ str(show_H)}.png")   # ergänzt um Charge

    # reference lines and legend panels
    fig = plt.gcf(); 
    fig.set_size_inches(10, 6)
    if config.name == "E333": fig.set_size_inches(20, 12)        
    fig.savefig(output_png, dpi=200)
    open_image(output_png)

    #     Legende
    add_legend_panels(Obj, i_Emax, D_i_c_)
    fig = plt.gcf(); 
    output_png = os.path.join(RESULTS_DIR, f"{config.name +"_legend" }.png")   
    fig.savefig(output_png, dpi=100)
    open_image(output_png)

    end_time_stamp = datetime.datetime.now()
    delta = end_time_stamp - start_time_stamp
    print(f"took {delta.seconds} seconds")

if __name__ == "__main__":

    main()

