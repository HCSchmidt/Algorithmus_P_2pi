from __future__ import annotations

import datetime
from pathlib import Path

import matplotlib.pyplot as plt

from .cli import parse_args, ScanSector, select_preset_by_sector
from .output.config import RESULTS_DIR
from .particles import PARTICLES
from .energy import EnergieEngine
from .output.plotting import (
    draw_points,
    add_reference_lines,
    add_legend_panels,
    add_particle_labels,
)
from .output.report import write_report_csv
from .utils import open_image


def _legend_stats_from_results(results, Cnt):
    """
    Build i_Emax (total Δi) and D_i_c_ (Δi/(2π) %) dictionaries for the legend,
    matching the semantics of the original script.
    """
    i_Emax = {}
    D_i_c_ = {}

    for key, pdata in results.items():
        total_delta_i = 0
        total_counts = 0

        for m in pdata.m_values():
            if pdata.i_Emax[m] == 0:
                continue

            delta_i = pdata.i_Emax[m] - pdata.i_Emin[m] + 1
            total_delta_i += abs(delta_i)

            c = Cnt[m]
            if c:
                total_counts += c

        i_Emax[key] = total_delta_i
        if total_counts > 0:
            D_i_c_[key] = f"{(total_delta_i * 100 / total_counts):.5f} %"
        else:
            D_i_c_[key] = ""

    return i_Emax, D_i_c_


def main(argv=None):
    start_time = datetime.datetime.now()

    args = parse_args(argv)
    sector: ScanSector = args.sector
    no_show: bool = args.no_show

    config = select_preset_by_sector(sector)

    results_dir = Path(RESULTS_DIR)
    results_dir.mkdir(parents=True, exist_ok=True)

    base_name = config.name
    png_path = results_dir / f"{base_name}.png"
    csv_path = results_dir / f"{base_name}.csv"

    engine = EnergieEngine(sector=sector)

    (
        xs_by_particle,
        ys_by_particle,
        grey_segments,
        results,
        Cnt,
        labels,
        i_T,
        i_T1,
    ) = engine.run(PARTICLES)  # OK: engine.run() now accepts dict

    print("possible ET:", i_T, "real ET:", i_T1)

    i_Emax, D_i_c_ = _legend_stats_from_results(results, Cnt)

    plt.figure(figsize=(10, 6))

    draw_points(xs_by_particle, ys_by_particle, grey_segments)
    add_reference_lines(sector, i_T)
    add_legend_panels(sector, i_T, PARTICLES, i_Emax, D_i_c_)
    add_particle_labels(labels)

    plt.tight_layout()
    plt.savefig(png_path, dpi=120)
    plt.close()


    write_report_csv(csv_path, sector, PARTICLES, results, Cnt)

    if not no_show:
        open_image(str(png_path))

    end_time = datetime.datetime.now()
    print(f"Finished in {(end_time - start_time).seconds} seconds")
    print(f"PNG  → {png_path}")
    print(f"CSV  → {csv_path}")


if __name__ == "__main__":
    main()