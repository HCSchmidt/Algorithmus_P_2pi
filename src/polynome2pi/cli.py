from enum import Enum
from dataclasses import dataclass
import argparse

class ScanSector(str, Enum):
    minimal = "001 u d"
    light = "011 u d s"
    broad = "112 u d s Nukleon"
    nucleon = "111 H-Atom"
    heavy = "222 u d s c"
    E112P = "112 E 1700

@dataclass
class PolynomeConfig:
    J4: int
    J3: int
    J2: int
    add_info: str = ""

    @property
    def name(self) -> str:
        suffix = f"_{self.add_info}" if self.add_info else ""
        return f"Polynom_{self.J4}{self.J3}{self.J2}{suffix}"


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description="Run P(2π) polynomial scan and plot results")

    parser.add_argument(
        "--sector",
        type=ScanSector,
        choices=list(ScanSector),
        required=True,
        help=(
            "Physical scan sector (polynomial depth). "
            f"Must be one of: {', '.join([e.value for e in ScanSector])}."
        ),
    )

    parser.add_argument(
        "--no-show",
        action="store_true",
        help="Do not open a GUI window; still save the PNG output",
    )


    return parser.parse_args(argv)


def select_preset_by_sector(sector: ScanSector) -> PolynomeConfig:
    if sector is ScanSector.minimal:
        return PolynomeConfig(J4=0, J3=0, J2=1, add_info=sector.value)

    if sector is ScanSector.light:
        return PolynomeConfig(J4=0, J3=1, J2=1, add_info=sector.value)

    if sector is ScanSector.broad:
        return PolynomeConfig(J4=1, J3=1, J2=2, add_info=sector.value)

    if sector is ScanSector.nucleon:
        return PolynomeConfig(J4=1, J3=1, J2=1, add_info=sector.value)

    if sector is ScanSector.heavy:
        return PolynomeConfig(J4=2, J3=2, J2=2, add_info=sector.value)

    if sector is ScanSector.E112P:
        return PolynomeConfig(J4=1, J3=1, J2=2, add_info=sector.value)    


    raise ValueError(f"Unhandled sector: {sector}")
