from enum import Enum
from dataclasses import dataclass
import argparse

class ScanSector(str, Enum):
    minimal = "000 u d"
    light = "012 u d s"
    broad = "112 u d s nucleon"
    nucleon = "111 H-atom"
    heavy = "222 u d s c"
    E112P = "112 E > 1700"
    E333U = "333 E > 1700"  
    E333D = "333 E < 2000" 
    E222 = "222 u d s c E < 3000" 

@dataclass
class PolynomeConfig:
    J4: int
    J3: int
    J2: int
    add_info: int

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
        return PolynomeConfig(J4=0, J3=0, J2=0, add_info=1)

    if sector is ScanSector.light:
        return PolynomeConfig(J4=0, J3=1, J2=2, add_info=2)

    if sector is ScanSector.broad:
        return PolynomeConfig(J4=1, J3=1, J2=2, add_info=3)   

    if sector is ScanSector.nucleon:
        return PolynomeConfig(J4=1, J3=1, J2=1, add_info=4)  

    if sector is ScanSector.E112P:
        return PolynomeConfig(J4=1, J3=1, J2=2, add_info=5) 
    
    if sector is ScanSector.heavy:
        return PolynomeConfig(J4=2, J3=2, J2=2, add_info=6)   

    if sector is ScanSector.E333U:
        return PolynomeConfig(J4=3, J3=3, J2=3, add_info=7)    

    if sector is ScanSector.E333D:
        return PolynomeConfig(J4=3, J3=3, J2=3, add_info=8) 
    
    if sector is ScanSector.E222:
        return PolynomeConfig(J4=2, J3=2, J2=2, add_info=9) 

    raise ValueError(f"Unhandled sector: {sector}")
