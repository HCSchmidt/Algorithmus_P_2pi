from __future__ import annotations

from dataclasses import dataclass
from .cli import ScanSector


@dataclass(frozen=True)
class ScanPreset:
    sector: ScanSector
    J4: int
    J3: int
    J2: int
    name: str


def preset_for_sector(sector: ScanSector) -> ScanPreset:
    if sector is ScanSector.minimal:
        return ScanPreset(sector, J4=0, J3=0, J2=1, name="minimal")
    if sector is ScanSector.light:
        return ScanPreset(sector, J4=0, J3=1, J2=1, name="light")
    if sector is ScanSector.broad:
        return ScanPreset(sector, J4=1, J3=1, J2=2, name="broad")
    if sector is ScanSector.nucleon:
        return ScanPreset(sector, J4=1, J3=1, J2=1, name="nucleon")
    if sector is ScanSector.heavy:
        return ScanPreset(sector, J4=2, J3=2, J2=2, name="heavy")
    raise ValueError(f"Unhandled sector: {sector}")
