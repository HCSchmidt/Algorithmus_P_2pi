from __future__ import annotations

from pathlib import Path
import os
from .cli import parse_args, ScanSector
from .presets import preset_for_sector
from .particles import get_particles  # dict[str, Particle]
from .utils import open_image
from .engine import ScanEngine
from .energy_model import EnergyModel
from .report import write_results_csv
from .plotting import plot_scan

from .sensitivity import run_sensitivity, write_sensitivity_csv, plot_sensitivity


def main(argv=None) -> int:
    args = parse_args(argv)
    sector: ScanSector = args.sector

    particles = get_particles()  # <-- wichtig: einmal holen, überall verwenden



    
    results_dir = Path(os.path.join("results", sector.value)).resolve()
    results_dir.mkdir(parents=True, exist_ok=True)

    preset = preset_for_sector(sector)

    # ---------------------------
    # Sensitivity mode
    # ---------------------------
    if args.command == "sensitivity":
        eps = args.eps
        steps = args.steps
        particle_key = args.particle
        hit_mode = args.hit_mode
        if particle_key is not None and particle_key not in particles:
            available = ", ".join(sorted(particles.keys()))
            raise SystemExit(f"Unknown particle key '{particle_key}'. Available keys: {available}")

        points = run_sensitivity(
            preset=preset,
            particles=particles,
            particle_key=particle_key,
            eps=float(eps),
            steps=int(steps),
            hit_mode=hit_mode,
        )

        suffix = f"_{particle_key}" if particle_key else ""
        csv_path = results_dir / f"sensitivity_{eps}_{sector.value}{suffix}.csv"
        png_path = results_dir / f"sensitivity_{eps}_{sector.value}{suffix}.png"

        # ✅ richtig: Sensitivity-CSV
        write_sensitivity_csv(points, csv_path)

        if particle_key:
            y_left = "particle_hits"
            y_right = "particle_hit_ratio"
            title = f"Sensitivity – sector={sector.value}, particle={particle_key}"
        else:
            y_left = "real_ET"
            y_right = "hit_ratio"
            title = f"Sensitivity – sector={sector.value} (global)"

        plot_sensitivity(
            points,
            png_path,
            title=title,
            y_left=y_left,
            y_right=y_right,
        )

        if not args.no_open:
            open_image(str(png_path))

        return 0

    # ---------------------------
    # Normal scan mode
    # ---------------------------
    model = EnergyModel()
    engine = ScanEngine(preset=preset, model=model)

    outputs = engine.run(particles)

    base_name = f"scan_{sector.value}"
    png_path = results_dir / f"{base_name}.png"
    csv_path = results_dir / f"{base_name}.csv"
    title=f"P(2π) scan – sector: {preset.name}",
    plot_scan(out_png=png_path, particles=particles, matched_points=outputs.matched_points,
              unmatched_segments=outputs.unmatched_segments,
              title=title)

    # ✅ richtig: write_results_csv erwartet bins_by_particle, nicht outputs
    write_results_csv(
        path=csv_path,
        sector=sector,
        particles=particles,
        bins_by_particle=outputs.bins_by_particle,
    )

    if not args.no_open:
        open_image(str(png_path))

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
