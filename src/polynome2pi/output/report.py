from __future__ import annotations

import csv
from pathlib import Path
from typing import Dict, Iterable, List

from ..cli import ScanSector
from ..particles import Particle


CSV_HEADER = [
    "sector",
    "particle_key",
    "particle_name",
    "symbol",
    "block",
    "row_type",
    "theory_E",
    "E_value",
    "i_value",
    "i4",
    "i3",
    "i2",
    "i1",
    "i0",
    "i_minus1",
    "counts",
    "delta_i",
    "delta_i_over_2pi",
]


def write_report_csv(
    file_path: Path,
    sector: ScanSector,
    particles: Dict[str, Particle],
    results,
    Cnt,
):
    """
    Write a single CSV containing all report data.

    `results` is expected to be a dict keyed by particle.key with:
        - Emax[m], Emin[m]
        - i_Emax[m], i_Emin[m]
        - Dmax[m, k], Dmin[m, k]
    """

    rows: List[List] = []

    for particle_key, particle in particles.items():
        pdata = results.get(particle_key)
        if pdata is None:
            continue

        block_index = 0
        total_delta_i = 0
        total_counts = 0

        for m in pdata.m_values():
            if pdata.i_Emax[m] == 0:
                continue

            block_index += 1

            E_max = pdata.Emax[m]
            E_min = pdata.Emin[m]
            E_mean = (E_max + E_min) / 2

            delta_i = pdata.i_Emax[m] - pdata.i_Emin[m] + 1
            counts = Cnt[m]
            delta_i_over_2pi = abs(delta_i) * 100 / counts

            total_delta_i += abs(delta_i)
            total_counts += counts

            # -------------------- max --------------------
            rows.append(_row(
                sector, particle, block_index, "max",
                particle.theory_E, E_max, pdata.i_Emax[m],
                pdata.Dmax[m], counts, delta_i, delta_i_over_2pi
            ))

            # -------------------- mean -------------------
            rows.append(_row(
                sector, particle, block_index, "mean",
                particle.theory_E, E_mean, delta_i,
                None, counts, delta_i, delta_i_over_2pi
            ))

            # -------------------- min --------------------
            rows.append(_row(
                sector, particle, block_index, "min",
                particle.theory_E, E_min, pdata.i_Emin[m],
                pdata.Dmin[m], counts, delta_i, delta_i_over_2pi
            ))

            # -------------------- delta ------------------
            rows.append(_row(
                sector, particle, block_index, "delta",
                particle.theory_E, "",
                delta_i, None, counts, delta_i, delta_i_over_2pi
            ))

        # -------------------- total --------------------
        if total_counts > 0:
            rows.append([
                sector.value,
                particle.key,
                particle.name,
                particle.symbol,
                "",
                "total",
                particle.theory_E,
                "",
                total_delta_i,
                "", "", "", "", "", "",
                total_counts,
                total_delta_i,
                total_delta_i * 100 / total_counts,
            ])

        # -------------------- info (only with i4 > 1) ---
        if block_index == 0:
            rows.append([
                sector.value,
                particle.key,
                particle.name,
                particle.symbol,
                "",
                "info",
                particle.theory_E,
                "",
                "",
                "", "", "", "", "", "",
                "",
                "",
                "only with i4 > 1",
            ])

    _write_csv(file_path, rows)


# ---------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------

def _row(
    sector: ScanSector,
    particle: Particle,
    block: int,
    row_type: str,
    theory_E,
    E_value,
    i_value,
    D_vals,
    counts,
    delta_i,
    delta_i_over_2pi,
):
    if D_vals is None:
        D_vals = ["", "", "", "", "", ""]

    return [
        sector.value,
        particle.key,
        particle.name,
        particle.symbol,
        block,
        row_type,
        theory_E,
        E_value,
        i_value,
        D_vals[0] / 2 if D_vals[0] != "" else "",
        D_vals[1] / 2 if D_vals[1] != "" else "",
        D_vals[2] / 2 if D_vals[2] != "" else "",
        D_vals[3] / 2 if D_vals[3] != "" else "",
        D_vals[4] / 2 if D_vals[4] != "" else "",
        D_vals[5] / 2 if D_vals[5] != "" else "",
        counts,
        delta_i,
        round(delta_i_over_2pi, 6),
    ]


def _write_csv(path: Path, rows: Iterable[List]):
    path.parent.mkdir(parents=True, exist_ok=True)

    with path.open("w", newline="", encoding="utf8") as f:
        writer = csv.writer(f)
        writer.writerow(CSV_HEADER)
        writer.writerows(rows)