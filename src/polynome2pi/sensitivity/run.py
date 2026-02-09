from typing import Optional
import numpy as np
from typing import Dict, List

from polynome2pi.cli import ScanSector
from polynome2pi.energy_model import EnergyModel
from polynome2pi.engine import ScanEngine, ScanOutputs
from polynome2pi.particles import Particle, get_particles
from polynome2pi.presets import ScanPreset, preset_for_sector
from polynome2pi.sensitivity.SensitivityPoint import SensitivityPoint
from polynome2pi.sensitivity.plotting import plot_sensitivity
from polynome2pi.sensitivity.report import write_sensitivity_csv


from pathlib import Path


def _particle_hits_from_outputs(outputs: ScanOutputs, particle_key: str, mode: str) -> int:
    """
    mode:
      - "matched_points": number of matched points for the particle (len(xs))
      - "delta_sum": sum of delta_i over all bins for the particle (more "legacy-like")
    """
    if mode == "matched_points":
        xs, _ys = outputs.matched_points.get(particle_key, ([], []))
        return int(len(xs))

    if mode == "delta_sum":
        bins = outputs.bins_by_particle.get(particle_key, [])
        return int(sum(b.delta_i for b in bins))

    raise ValueError(f"Unknown particle hit mode: {mode}")


def run_sensitivity(
    preset: ScanPreset,
    particles: Dict[str, Particle],
    *,
    particle_key: Optional[str],
    eps: float,
    steps: int,
    hit_mode: str = "matched_points",
) -> List[SensitivityPoint]:
    """
    Sweep base_scale around 1.0 (== exact 2π), run full scan each time,
    and collect global + optional per-particle metrics.
    """
    scales = np.linspace(1.0 - eps, 1.0 + eps, int(steps))

    points: List[SensitivityPoint] = []
    for s in scales:
        model = EnergyModel(base_scale=float(s))  # <--- benötigt EnergyModel(base_scale=...)
        engine = ScanEngine(preset=preset, model=model)

        outputs = engine.run(particles)

        accepted_scan_points = int(outputs.accepted_scan_points) if outputs.accepted_scan_points else 0
        total_particle_hits = int(outputs.total_particle_hits)
        hit_ratio = (total_particle_hits / accepted_scan_points) if accepted_scan_points else 0.0

        if particle_key:
            phits = _particle_hits_from_outputs(outputs, particle_key, mode=hit_mode)
            phit_ratio = (phits / accepted_scan_points) if accepted_scan_points else 0.0
            points.append(
                SensitivityPoint(
                    base_scale=float(s),
                    accepted_scan_points=accepted_scan_points,
                    total_particle_hits=total_particle_hits,
                    hit_ratio=hit_ratio,
                    particle_key=particle_key,
                    particle_hits=phits,
                    particle_hit_ratio=phit_ratio,
                )
            )
        else:
            points.append(
                SensitivityPoint(
                    base_scale=float(s),
                    accepted_scan_points=accepted_scan_points,
                    total_particle_hits=total_particle_hits,
                    hit_ratio=hit_ratio,
                )
            )

    return points


def evaluate_sensitivity(
    sector: ScanSector, eps: float, steps: int, particle_key: str, hit_mode: str, results_dir: Path
) -> Path:
    particles = get_particles()
    if particle_key is not None and particle_key not in particles:
        available = ", ".join(sorted(particles.keys()))
        raise SystemExit(f"Unknown particle key '{particle_key}'. Available keys: {available}")

    preset = preset_for_sector(sector)
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

    write_sensitivity_csv(points, csv_path)

    if particle_key:
        y_left = "particle_hits"
        y_right = "particle_hit_ratio"
        title = f"Sensitivity – sector={sector.value}, particle={particle_key}"
    else:
        y_left = "total_particle_hits"
        y_right = "hit_ratio"
        title = f"Sensitivity – sector={sector.value} (global)"

    plot_sensitivity(
        points,
        png_path,
        title=title,
        y_left=y_left,
        y_right=y_right,
    )

    return png_path
