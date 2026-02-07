import argparse
from enum import Enum


class ScanSector(str, Enum):
    minimal = "minimal"
    light = "light"
    broad = "broad"
    nucleon = "nucleon"
    heavy = "heavy"


def parse_args(argv=None):
    # Parent parser for shared args
    parent = argparse.ArgumentParser(add_help=False)
    parent.add_argument(
        "--sector",
        type=ScanSector,
        choices=list(ScanSector),
        required=True,
        help="Scan sector (preset).",
    )
    parent.add_argument(
        "--no-open",
        action="store_true",
        help="Do not open the generated PNG automatically.",
    )

    parser = argparse.ArgumentParser(
        description="Scan P(2π) coefficient space and match particle energies."
    )

    sub = parser.add_subparsers(
        dest="command",
        required=True,
        metavar="{scan,sensitivity}",
    )

    # -----------------
    # scan command
    # -----------------
    scan_p = sub.add_parser(
        "scan",
        parents=[parent],
        help="Run a normal scan and produce plot + CSV.",
    )
    # (Optional: scan-specific args could go here)

    # -----------------
    # sensitivity command
    # -----------------
    sens_p = sub.add_parser(
        "sensitivity",
        parents=[parent],
        help="Run sensitivity analysis by varying the base around 2π.",
    )

    sens_p.add_argument(
        "--eps",
        type=float,
        default=0.001,
        help="Relative variation around 2π (±eps). Default: 0.001 (±0.1%%).",
    )

    sens_p.add_argument(
        "--steps",
        type=int,
        default=21,
        help="Number of base samples in the sensitivity sweep. Default: 21.",
    )

    sens_p.add_argument(
        "--particle",
        type=str,
        default=None,
        help="If set, run sensitivity for a single particle key (e.g. 'muon', 'proton', 'pion_0', ...).",
    )

    sens_p.add_argument(
        "--hit-mode",
        type=str,
        default="matched_points",
        choices=["matched_points", "delta_sum"],
        help=(
            "How to count particle hits for sensitivity. "
            "'matched_points' = number of matched scan points; "
            "'delta_sum' = sum of delta_i over bins (more legacy-like)."
        ),
    )

    return parser.parse_args(argv)