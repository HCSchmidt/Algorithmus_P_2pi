from __future__ import annotations

from pathlib import Path
import os

from polynome2pi.sensitivity.run import evaluate_sensitivity
from polynome2pi.scan.run import run_scan
from .cli import parse_args, ScanSector, EvalulationType
from .utils import open_image


def main(argv=None) -> int:
    args = parse_args(argv)
    sector: ScanSector = args.sector

    results_dir = Path(os.path.join("results", sector.value)).resolve()
    results_dir.mkdir(parents=True, exist_ok=True)

    if args.command == EvalulationType.SENSITIVITY:

        png_path = evaluate_sensitivity(
            sector=sector,
            eps=args.eps,
            steps=args.steps,
            particle_key=args.particle,
            hit_mode=args.hit_mode,
            results_dir=results_dir,
        )
        if not args.no_open:
            open_image(png_path)

        return 0

    if args.command == EvalulationType.SCAN:

        png_path, png_grid_path = run_scan(sector, results_dir)

        if not args.no_open:
            open_image(png_path)
            open_image(png_grid_path)

        return 0


if __name__ == "__main__":
    raise SystemExit(main())
