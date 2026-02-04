from __future__ import annotations

import datetime
import subprocess
import sys
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


    # ------------------------------------------------------------------
    # Energy scan
    # ------------------------------------------------------------------
    engine = EnergieEngine(sector=sector)

    (
        xs_by_j,
        ys_by_j,
        grey_segments,
        results,
        Cnt,
        labels,
        i_T,
        i_T1,
    ) = engine.run(PARTICLES)

    print("possible ET:", i_T, "real ET:", i_T1)

    # ------------------------------------------------------------------
    # Plot
    # ------------------------------------------------------------------
    plt.figure(figsize=(10, 6))

    draw_points(xs_by_j, ys_by_j, grey_segments)
    add_reference_lines(sector, i_T)
    add_legend_panels(sector, i_T, PARTICLES, results)
    add_particle_labels(labels)

    plt.tight_layout()
    plt.savefig(png_path, dpi=120)
    plt.close()

    # ------------------------------------------------------------------
    # CSV report
    # ------------------------------------------------------------------
    write_report_csv(
        file_path=csv_path,
        sector=sector,
        particles=PARTICLES,
        results=results,
        Cnt=Cnt,
    )

    # ------------------------------------------------------------------
    # Open PNG (OS-independent)
    # ------------------------------------------------------------------
    if not no_show:
        open_file_cross_platform(png_path)

    end_time = datetime.datetime.now()
    print(f"Finished in {(end_time - start_time).seconds} seconds")
    print(f"PNG  → {png_path}")
    print(f"CSV  → {csv_path}")


if __name__ == "__main__":
    main()