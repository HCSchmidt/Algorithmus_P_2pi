from dataclasses import dataclass
from typing import Optional


@dataclass(frozen=True)
class SensitivityPoint:
    base_scale: float
    accepted_scan_points: int
    total_particle_hits: int
    hit_ratio: float

    # optional per-particle metrics
    particle_key: Optional[str] = None
    particle_hits: Optional[int] = None  # count of matched points (or delta sum, siehe unten)
    particle_hit_ratio: Optional[float] = None  # particle_hits / possible_ET
