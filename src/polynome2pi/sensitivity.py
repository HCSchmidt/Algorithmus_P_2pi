from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import matplotlib.pyplot as plt

from .engine import ScanEngine, ScanOutputs
from .energy_model import EnergyModel
from .particles import Particle
from .presets import ScanPreset


@dataclass(frozen=True)
class SensitivityPoint:
    base_scale: float
    possible_ET: int
    real_ET: int
    hit_ratio: float

    # optional per-particle metrics
    particle_key: Optional[str] = None
    particle_hits: Optional[int] = None  # count of matched points (or delta sum, siehe unten)
    particle_hit_ratio: Optional[float] = None  # particle_hits / possible_ET


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

        possible_ET = int(outputs.possible_ET) if outputs.possible_ET else 0
        real_ET = int(outputs.real_ET)
        hit_ratio = (real_ET / possible_ET) if possible_ET else 0.0

        if particle_key:
            phits = _particle_hits_from_outputs(outputs, particle_key, mode=hit_mode)
            phit_ratio = (phits / possible_ET) if possible_ET else 0.0
            points.append(
                SensitivityPoint(
                    base_scale=float(s),
                    possible_ET=possible_ET,
                    real_ET=real_ET,
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
                    possible_ET=possible_ET,
                    real_ET=real_ET,
                    hit_ratio=hit_ratio,
                )
            )

    return points


def write_sensitivity_csv(points: List[SensitivityPoint], csv_path: Path) -> None:
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with csv_path.open("w", encoding="utf-8") as f:
        # columns intentionally simple + stable
        f.write(
            "base_scale,possible_ET,real_ET,hit_ratio,particle_key,particle_hits,particle_hit_ratio\n"
        )
        for p in points:
            f.write(
                f"{p.base_scale:.8f},"
                f"{p.possible_ET},"
                f"{p.real_ET},"
                f"{p.hit_ratio:.10f},"
                f"{p.particle_key or ''},"
                f"{'' if p.particle_hits is None else p.particle_hits},"
                f"{'' if p.particle_hit_ratio is None else f'{p.particle_hit_ratio:.10f}'}\n"
            )


def plot_sensitivity(
    points: List[SensitivityPoint],
    png_path: Path,
    *,
    title: str,
    y_left: str = "real_ET",
    y_right: str = "hit_ratio",
) -> None:
    """
    y_left:
      - "real_ET" (global)
      - "particle_hits" (only if particle_key set)
    """
    xs = [p.base_scale for p in points]

    if y_left == "real_ET":
        y1 = [p.real_ET for p in points]
        y1_label = "real_ET"
    elif y_left == "particle_hits":
        y1 = [p.particle_hits or 0 for p in points]
        y1_label = "particle_hits"
    else:
        raise ValueError(f"Unknown y_left: {y_left}")

    if y_right == "hit_ratio":
        y2 = [p.hit_ratio for p in points]
        y2_label = "hit_ratio"
    elif y_right == "particle_hit_ratio":
        y2 = [p.particle_hit_ratio or 0.0 for p in points]
        y2_label = "particle_hit_ratio"
    else:
        raise ValueError(f"Unknown y_right: {y_right}")

    png_path.parent.mkdir(parents=True, exist_ok=True)

    fig, ax1 = plt.subplots(figsize=(10, 6))
    ax1.plot(xs, y1, marker="o")
    ax1.set_xlabel("Base scale (relative to 2π)")
    ax1.set_ylabel(y1_label)
    ax1.grid(True)

    ax2 = ax1.twinx()
    ax2.plot(xs, y2, linestyle="--")
    ax2.set_ylabel(y2_label)

    plt.title(title)
    plt.tight_layout()
    plt.savefig(png_path, dpi=120)
    plt.close(fig)
