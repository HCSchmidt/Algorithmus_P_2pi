import os
import csv
from pathlib import Path

from ..cli import ScanSector
from ..constants import get_colors
from .config import RESULTS_DIR

def write_report(
    base_name: str,
    sector: ScanSector,
    Obj,
    D_i_N,
    D_i_c_,
    Cnt,
    Emax,
    Emin,
    i_Emax,
    i_Emin,
    Dmax,
    Dmin,
):
    """
    Writes two CSV files:
      - <stem>_events.csv  : one row per (particle, m) hit, containing max/mean/min + deltas
      - <stem>_summary.csv : one row per particle, containing total ΣΔi, total counts, percent

    Returns:
      labels: list of tuples (sector, j, x_base, y, text, color) used for plotting.
    """
    colors = get_colors()
    labels = []

    events_rows = []
    summary_rows = []

    sector_value = getattr(sector, "value", str(sector))

    for j in range(1, 29):
        # theory bounds (same as TXT)
        m_min = float(Obj[j][2]) + float(Obj[j][3])
        m_max = float(Obj[j][2]) + float(Obj[j][4])
        p_g = len(Obj[j][2])
        m_min = float(str(m_min)[:p_g])
        m_max = float(str(m_max)[:p_g])

        had_hits = False
        labeled_this_particle = False

        for m in range(1, 512):
            if i_Emax[j, m] == 0:
                continue

            had_hits = True

            # same calculations as TXT
            E_mean = (Emax[j, m] + Emin[j, m]) / 2
            Di_E = i_Emax[j, m] - i_Emin[j, m] + 1
            i_Emax[j, 0] += abs(Di_E)

            Cnt_ = Cnt[m]
            D_i_c = round((abs(Di_E)) * 100 / Cnt_, 5)
            i_Emax[j, 516] += Cnt_
            D_i_N[j] = float(i_Emax[j, 0]) * 100 / Cnt[m]

            # truncation like TXT
            Emax_val = float(str(Emax[j, m])[:p_g])
            Emin_val = float(str(Emin[j, m])[:p_g])
            E_mean_val = float(str(E_mean)[:p_g])

            # plot label once per particle (first hit)
            if not labeled_this_particle:
                labels.append((sector, j, i_Emin[j, m], E_mean_val, Obj[j][9], colors[j]))
                labeled_this_particle = True

            # ---- events CSV row (wide format: max/mean/min + deltas + coeffs for max & min) ----
            events_rows.append({
                "sector": sector_value,
                "particle": Obj[j][0],
                "particle_symbol": Obj[j][9],
                "m_bin": m,

                "theory_E": Obj[j][2],
                "theory_m_min": m_min,
                "theory_m_max": m_max,

                "E_max": Emax_val,
                "E_mean": E_mean_val,
                "E_min": Emin_val,

                "N_max": int(i_Emax[j, m]),
                "N_min": int(i_Emin[j, m]),

                "delta_i": int(Di_E),
                "abs_delta_i": int(abs(Di_E)),
                "counts": int(Cnt_),
                "delta_i_percent": float(D_i_c),

                # Dmax coefficients
                "max_i4": Dmax[j, m, 0] / 2,
                "max_i3": Dmax[j, m, 1] / 2,
                "max_i2": Dmax[j, m, 2] / 2,
                "max_i1": Dmax[j, m, 3] / 2,
                "max_i0": Dmax[j, m, 4] / 2,
                "max_i_minus_1": Dmax[j, m, 5] / 2,
                "max_C": Dmax[j, m, 6] / 2,

                # Dmin coefficients
                "min_i4": Dmin[j, m, 0] / 2,
                "min_i3": Dmin[j, m, 1] / 2,
                "min_i2": Dmin[j, m, 2] / 2,
                "min_i1": Dmin[j, m, 3] / 2,
                "min_i0": Dmin[j, m, 4] / 2,
                "min_i_minus_1": Dmin[j, m, 5] / 2,
                "min_C": Dmin[j, m, 6] / 2,
            })

        # ---- summary row per particle (like TXT "total" line) ----
        if j > 1 and had_hits and i_Emax[j, 516] > 0:
            D_i_c_tot = round(float(i_Emax[j, 0]) * 100 / i_Emax[j, 516], 5)
            D_i_c_[j] = f"{D_i_c_tot} %"

            summary_rows.append({
                "sector": sector_value,
                "particle": Obj[j][0],
                "particle_symbol": Obj[j][9],
                "theory_E": Obj[j][2],
                "sum_abs_delta_i": int(i_Emax[j, 0]),
                "total_counts": int(i_Emax[j, 516]),
                "delta_i_percent": float(D_i_c_tot),
            })

    # ---- write events CSV ----
    if events_rows:
        events_path = os.path.join(RESULTS_DIR, f"{base_name}_events.csv")
        with open(events_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=list(events_rows[0].keys()))
            writer.writeheader()
            writer.writerows(events_rows)

    # ---- write summary CSV ----
    if summary_rows:
        summary_path = os.path.join(RESULTS_DIR, f"{base_name}_summary.csv")
        with open(summary_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=list(summary_rows[0].keys()))
            writer.writeheader()
            writer.writerows(summary_rows)

    return labels
