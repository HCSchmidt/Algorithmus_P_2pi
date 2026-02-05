from __future__ import annotations

import argparse
from enum import Enum


class ScanSector(str, Enum):
    minimal = "minimal"
    light = "light"
    broad = "broad"
    nucleon = "nucleon"
    heavy = "heavy"


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description="Scan P(2π) coefficient space and match particle energies.")

    parser.add_argument(
        "--sector",
        type=ScanSector,
        choices=list(ScanSector),
        required=True,
        help="Scan preset. One of: minimal/light/broad/nucleon/heavy",
    )

    parser.add_argument(
        "--out-dir",
        default="results",
        help="Output directory for CSV and PNG (default: results)",
    )

    parser.add_argument(
        "--no-open",
        action="store_true",
        help="Do not open the generated PNG automatically.",
    )

    return parser.parse_args(argv)