from __future__ import annotations

from pathlib import Path
import time

from .cli import parse_args
from .presets import preset_for_sector
from .particles import get_particles
from .energy_model import EnergyModel
from .engine import ScanEngine
from .report import write_results_csv
from .plotting import plot_scan, open_file
from .sensitivity import run_sensitivity, write_sensitivity_csv, plot_sensitivity

def main(argv=None) -> int:
    args = parse_args(argv)
    preset = preset_for_sector(args.sector)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    
    particles = get_particles()

    if args.sensitivity:
        rows = run_sensitivity(
            sector=args.sector,
            preset=preset,
            particles=particles,
            eps=args.eps,
            steps=args.steps,
        )
        csv_path = out_dir / f"results_{preset.name}_sensitivity_eps_{args.eps}.csv"
        png_path = out_dir / f"plot_{preset.name}_sensitivity_eps_{args.eps}.png"

        write_sensitivity_csv(csv_path, rows)
        plot_sensitivity(rows, png_path)
        return

    model = EnergyModel()
    engine = ScanEngine(preset, model)

    outputs = engine.run(particles)

    print(f"possible ET: {outputs.possible_ET}  real ET: {outputs.real_ET}")

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