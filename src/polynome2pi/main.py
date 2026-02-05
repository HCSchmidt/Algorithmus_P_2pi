from __future__ import annotations

from pathlib import Path
import time

from .cli import parse_args
from .presets import preset_for_sector
from .particles import default_particles
from .energy_model import EnergyModel
from .engine import ScanEngine
from .report import write_results_csv
from .plotting import plot_scan, open_file


def main(argv=None) -> int:
    args = parse_args(argv)
    preset = preset_for_sector(args.sector)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    particles = default_particles()
    model = EnergyModel()
    engine = ScanEngine(preset, model)

    t0 = time.time()
    outputs = engine.run(particles)
    dt = time.time() - t0

    print(f"possible ET: {outputs.possible_ET}  real ET: {outputs.real_ET}")
    print(f"took {dt:.2f} seconds")

    csv_path = out_dir / f"results_{preset.name}.csv"
    png_path = out_dir / f"plot_{preset.name}.png"

    write_results_csv(csv_path, preset.sector, particles, outputs.bins_by_particle)
    plot_scan(
        png_path,
        particles,
        outputs.matched_points,
        outputs.unmatched_segments,
        title=f"P(2π) scan – sector: {preset.name}",
    )

    print(f"Wrote CSV: {csv_path}")
    print(f"Wrote PNG: {png_path}")

    if not args.no_open:
        open_file(png_path)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())