from __future__ import annotations

import csv
from pathlib import Path
from typing import Dict, List

from polynome2pi.cli import ScanSector
from polynome2pi.particles import Particle
from polynome2pi.engine import BinResult


def write_results_csv(
    path: Path,
    sector: ScanSector,
    particles: Dict[str, Particle],
    bins_by_particle: Dict[str, List[BinResult]],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)

    headers = [
        "sector",
        "particle_key",
        "particle_name",
        "symbol",
        "theory_str",
        "theory_E",
        "sd_minus",
        "sd_plus",
        "m_bin",
        "E_mean",
        "E_min",
        "E_max",
        "i_T_min",
        "i_T_max",
        "delta_i",
        "counts",
        "delta_i_over_counts",
        # coefficients (stored raw ints; divide by 2 when interpreting)
        "min_i4",
        "min_i3",
        "min_i2",
        "min_i1",
        "min_i0",
        "min_i_minus_1",
        "min_C",
        "max_i4",
        "max_i3",
        "max_i2",
        "max_i1",
        "max_i0",
        "max_i_minus_1",
        "max_C",
    ]

    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter=";")
        w.writerow(headers)

        for pkey, p in particles.items():
            rows = bins_by_particle.get(pkey, [])
            if not rows:
                # still emit a summary row with no bins (helps downstream)
                w.writerow(
                    [
                        sector.value,
                        p.key,
                        p.name,
                        p.symbol,
                        p.theory_str,
                        p.theory_E,
                        p.sd_minus,
                        p.sd_plus,
                        "",
                        "",
                        "",
                        "",
                        "",
                        "",
                        "",
                        "",
                        "",  # bin fields
                        "",
                        "",
                        "",
                        "",
                        "",
                        "",
                        "",
                        "",
                        "",
                        "",
                        "",
                        "",
                        "",
                        "",
                        "",  # coeff fields
                    ]
                )
                continue

            for r in rows:
                (min_i4, min_i3, min_i2, min_i1, min_i0, min_im1, min_C) = r.coeff_min
                (max_i4, max_i3, max_i2, max_i1, max_i0, max_im1, max_C) = r.coeff_max

                w.writerow(
                    [
                        sector.value,
                        p.key,
                        p.name,
                        p.symbol,
                        p.theory_str,
                        p.theory_E,
                        p.sd_minus,
                        p.sd_plus,
                        r.m,
                        r.mean_E,
                        r.E_min,
                        r.E_max,
                        r.i_T_min,
                        r.i_T_max,
                        r.delta_i,
                        r.counts,
                        r.delta_i_over_counts,
                        min_i4,
                        min_i3,
                        min_i2,
                        min_i1,
                        min_i0,
                        min_im1,
                        min_C,
                        max_i4,
                        max_i3,
                        max_i2,
                        max_i1,
                        max_i0,
                        max_im1,
                        max_C,
                    ]
                )
