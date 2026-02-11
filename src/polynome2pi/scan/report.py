from __future__ import annotations

import csv
from pathlib import Path
from typing import Dict, List
import pandas as pd

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


def write_results_txt(
    path: Path,
    sector: ScanSector,
    particles: Dict[str, Particle],
    bins_by_particle: Dict[str, List[BinResult]],
) -> None:
    """
    Writes a human-readable legacy-like TXT report using pandas formatting (no manual spacing).

    Layout per particle:
      - particle name line
      - rows for each bin: max, mean, min, delta
      - optional total summary line
      - if no bins: a single 'mean ... only with i4 > 1' note line
    """
    path.parent.mkdir(parents=True, exist_ok=True)

    # ---------- format helpers ----------
    def fmt_energy(x):
        if x is None or (isinstance(x, float) and pd.isna(x)):
            return ""
        return f"{float(x):.11f}"

    def fmt_theory(x):
        # theory column sometimes holds strings like "1.0000...(31)"
        if x is None or (isinstance(x, float) and pd.isna(x)):
            return ""
        return str(x)

    def fmt_int(x):
        if x is None or (isinstance(x, float) and pd.isna(x)):
            return ""
        try:
            return f"{int(x)}"
        except Exception:
            return str(x)

    def fmt_coeff_half(x):
        if x is None or (isinstance(x, float) and pd.isna(x)):
            return ""
        return f"{(int(x) / 2.0):.1f}"

    def fmt_percent(x):
        if x is None or (isinstance(x, float) and pd.isna(x)):
            return ""
        return f"{float(x):.5f} %"

    header_cols = [
        "row_type",
        "theory",
        "E",
        "total",
        "i4",
        "i3",
        "i2",
        "i1",
        "i0",
        "i-1",
        "C",
        "cts",
        "di_over_2pi",
        "note",
    ]

    # ---------- write ----------
    with path.open("w", encoding="utf-8") as f:
        for pkey, p in particles.items():
            rows = bins_by_particle.get(pkey, [])

            # particle headline
            f.write(f"{p.name}\n")

            # no bins -> legacy note
            if not rows:
                # keep it close to old output:
                f.write(f"          mean  {p.theory_str}                    only with i4 > 1\n")
                continue

            table_rows: List[dict] = []

            sum_delta_i = 0
            sum_counts = 0

            for r in rows:
                (min_i4, min_i3, min_i2, min_i1, min_i0, min_im1, min_C) = r.coeff_min
                (max_i4, max_i3, max_i2, max_i1, max_i0, max_im1, max_C) = r.coeff_max

                theory_max = p.theory_E + p.sd_plus
                theory_min = p.theory_E + p.sd_minus

                # MAX
                table_rows.append(
                    dict(
                        row_type="max",
                        theory=theory_max,
                        E=r.E_max,
                        total=r.i_T_max,
                        i4=max_i4,
                        i3=max_i3,
                        i2=max_i2,
                        i1=max_i1,
                        i0=max_i0,
                        **{"i-1": max_im1},
                        C=max_C,
                        cts=None,
                        di_over_2pi=None,
                        note="",
                    )
                )

                # MEAN
                table_rows.append(
                    dict(
                        row_type="mean",
                        theory=p.theory_str,  # keep legacy theory string here
                        E=r.mean_E,
                        total=r.delta_i,  # legacy shows delta_i in this line under "total"
                        i4=None,
                        i3=None,
                        i2=None,
                        i1=None,
                        i0=None,
                        **{"i-1": None},
                        C=None,
                        cts=None,
                        di_over_2pi=None,
                        note="",
                    )
                )

                # MIN
                table_rows.append(
                    dict(
                        row_type="min",
                        theory=theory_min,
                        E=r.E_min,
                        total=r.i_T_min,
                        i4=min_i4,
                        i3=min_i3,
                        i2=min_i2,
                        i1=min_i1,
                        i0=min_i0,
                        **{"i-1": min_im1},
                        C=min_C,
                        cts=None,
                        di_over_2pi=None,
                        note="",
                    )
                )

                # DELTA (legacy-like: abs(delta_i), counts and percent)
                pct = 100.0 * r.delta_i_over_counts if r.counts else 0.0
                table_rows.append(
                    dict(
                        row_type="delta",
                        theory="",
                        E=abs(r.delta_i),
                        total="",  # keep clean
                        i4=None,
                        i3=None,
                        i2=None,
                        i1=None,
                        i0=None,
                        **{"i-1": None},
                        C=None,
                        cts=r.counts,
                        di_over_2pi=pct,
                        note="∆ abs(i)",
                    )
                )

                sum_delta_i += r.delta_i
                sum_counts += r.counts

            # TOTAL summary
            total_pct = 100.0 * (sum_delta_i / sum_counts) if sum_counts else 0.0
            table_rows.append(
                dict(
                    row_type="total",
                    theory="",
                    E=sum_delta_i,
                    total="Σ ∆i",
                    i4=None,
                    i3=None,
                    i2=None,
                    i1=None,
                    i0=None,
                    **{"i-1": None},
                    C=None,
                    cts=sum_counts,
                    di_over_2pi=total_pct,
                    note="",
                )
            )

            df = pd.DataFrame(table_rows, columns=header_cols)

            # apply formatting column-wise
            df["row_type"] = df["row_type"].map(lambda x: f"{x:>6}" if isinstance(x, str) else "")
            df["theory"] = df["theory"].map(fmt_theory)
            df["E"] = df["E"].map(fmt_energy)

            # "total" can be int or label
            df["total"] = df["total"].map(lambda x: fmt_int(x) if str(x).isdigit() else (x or ""))

            # coefficients: half + 1 decimal
            for c in ["i4", "i3", "i2", "i1", "i0", "i-1", "C"]:
                df[c] = df[c].map(fmt_coeff_half)

            df["cts"] = df["cts"].map(fmt_int)
            df["di_over_2pi"] = df["di_over_2pi"].map(fmt_percent)
            df["note"] = df["note"].fillna("")

            # render without index; pandas handles spacing/alignment
            txt = df.to_string(index=False, na_rep="")
            f.write(txt + "\n")
