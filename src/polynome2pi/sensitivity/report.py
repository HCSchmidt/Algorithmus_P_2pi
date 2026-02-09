from polynome2pi.sensitivity.SensitivityPoint import SensitivityPoint


from pathlib import Path
from typing import List


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
                f"{p.accepted_scan_points},"
                f"{p.total_particle_hits},"
                f"{p.hit_ratio:.10f},"
                f"{p.particle_key or ''},"
                f"{'' if p.particle_hits is None else p.particle_hits},"
                f"{'' if p.particle_hit_ratio is None else f'{p.particle_hit_ratio:.10f}'}\n"
            )
