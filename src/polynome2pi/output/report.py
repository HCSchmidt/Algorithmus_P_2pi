import os
import csv

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
    PEmax,
    PEmin,
):
    """
    Writes a single CSV file mirroring the old TXT table structure.   

    Returns:
      labels: list of tuples (sector, j, x_base, y, text, color) used for plotting.
    """
    colors = get_colors()
    labels = []

    rows = []

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

        Emax[j, 1028] = 0
        Emin[j, 1028] = 1000000
        for m in range(1, 1022):
            if i_Emax[j, m] != 0:
                had_hits = True
                break

        if not had_hits:
            rows.append({
                "particle": Obj[j][0],
                "row_type": "mean",
                "theory": Obj[j][1],
                "E": "",
                "total": "",
                "i4": "",
                "i3": "",
                "i2": "",
                "i1": "",
                "i0": "",
                "i_minus_1": "",
                "C": "",
                "cts": "",
                "di_over_2pi_percent": "",
                "note": "only with i4 > 1",
                "Polynom": "",
            })
            continue

        # Append particle header row
        rows.append({
            "particle": Obj[j][0],
            "row_type": "particle",
            "theory": "",
            "E": "",
            "total": "",
            "i4": "",
            "i3": "",
            "i2": "",
            "i1": "",
            "i0": "",
            "i_minus_1": "",
            "C": "",
            "cts": "",
            "di_over_2pi_percent": "",
            "note": "",
            "Polynom": ""
        })

        for m in range(1, 1022):
            if i_Emax[j, m] == 0:
                continue

            # same calculations as TXT
            E_mean = (Emax[j, m] + Emin[j, m]) / 2
            Di_E = i_Emax[j, m] - i_Emin[j, m] + 1
            i_Emax[j, 0] += abs(Di_E)

            Cnt_ = Cnt[m]
            D_i_c = round((abs(Di_E)) * 100 / Cnt_, 5)
            i_Emax[j, 1026] += Cnt_
            D_i_N[j] = float(i_Emax[j, 0]) * 100 / Cnt[m]

            # truncation like TXT
            Emax_val = float(str(Emax[j, m])[:p_g])
            Emin_val = float(str(Emin[j, m])[:p_g])
            E_mean_val = float(str(E_mean)[:p_g])

            if Emax[j, 1028] < Emax[j, m]:
                Emax[j, 1028] = Emax[j, m]; m__max = m
            if Emin[j, 1028] > Emin[j, m]:
                Emin[j, 1028] = Emin[j, m]; m__min = m

            # plot label once per particle (first hit)
            if not labeled_this_particle:
                labels.append((sector, j, i_Emin[j, m], E_mean_val, Obj[j][9], colors[j]))
                labeled_this_particle = True
 
            # max row
            rows.append({
                "particle": "",
                "row_type": "max",
                "theory": m_max,
                "E": Emax_val,
                "total": int(i_Emax[j, m]),
                "i4": Dmax[j, m, 0] / 2,
                "i3": Dmax[j, m, 1] / 2,
                "i2": Dmax[j, m, 2] / 2,
                "i1": Dmax[j, m, 3] / 2,
                "i0": Dmax[j, m, 4] / 2,
                "i_minus_1": Dmax[j, m, 5] / 2,
                "C": Dmax[j, m, 6] / 2,
                "cts": "",
                "di_over_2pi_percent": "",
                "note": "",
                "Polynom": PEmax[j][m] 
            })

            # mean row
            rows.append({
                "particle": "",
                "row_type": "mean",
                "theory": Obj[j][1],
                "E": E_mean_val,
                "total": int(Di_E),
                "i4": "",
                "i3": "",
                "i2": "",
                "i1": "",
                "i0": "",
                "i_minus_1": "",
                "C": "",
                "cts": "",
                "di_over_2pi_percent": "",
                "note": "",
                "Polynom": ""
            })

            # min row
            rows.append({
                "particle": "",
                "row_type": "min",
                "theory": m_min,
                "E": Emin_val,
                "total": int(i_Emin[j, m]),
                "i4": Dmin[j, m, 0] / 2,
                "i3": Dmin[j, m, 1] / 2,
                "i2": Dmin[j, m, 2] / 2,
                "i1": Dmin[j, m, 3] / 2,
                "i0": Dmin[j, m, 4] / 2,
                "i_minus_1": Dmin[j, m, 5] / 2,
                "C": Dmin[j, m, 6] / 2,
                "cts": "",
                "di_over_2pi_percent": "",
                "note": "",
                "Polynom": PEmin[j][m] 
            })

            # delta row
            rows.append({
                "particle": "",
                "row_type": "delta",
                "theory": "",
                "E": "",
                "total": int(abs(Di_E)),
                "i4": "",
                "i3": "",
                "i2": "",
                "i1": "",
                "i0": "",
                "i_minus_1": "",
                "C": "",
                "cts": int(Cnt_),
                "di_over_2pi_percent": float(D_i_c),
                "note": "∆ abs(i)",
                "Polynom": "",
            })

        # total row per particle if applicable
        if j > 1 and i_Emax[j, 1026] > 0:
            D_i_c_tot = round(float(i_Emax[j, 0]) * 100 / i_Emax[j, 1026], 5)
            D_i_c_[j] = f"{D_i_c_tot} %"

            for d in range(0,7):
                P= str(Dmin[j, m__max, 0]/2) + "(2 pi)^4 + "
                P += str(Dmin[j, m__max, 1]/2) + "(2 pi)^3 + "
                P += str(Dmin[j, m__max, 2]/2) + "(2 pi)^2 + "
                P += "("+ str(-Dmin[j, m__max, 3]/2) + " + " + str(-Dmin[j, m__max, 6]/2) +") (2 pi)^1 + "
                P += str(-Dmin[j, m__max, 4]/2) + " + "
                P += "(" + str(-Dmin[j, m__max, 5]/2) + " + " + str(Dmin[j, m__max, 6]) + ") (2 pi)^-1 "

            rows.append({
                "particle": "",
                "row_type": "total",
                "theory": "",
                "E": "",
                "E_m": Emax[j, 1028],
                "total": int(i_Emax[j, 0]),
                "i4": Dmin[j, m__max, 0] / 2,
                "i3": Dmin[j, m__max, 1] / 2,
                "i2": Dmin[j, m__max, 2] / 2,
                "i1": Dmin[j, m__max, 3] / 2,
                "i0": Dmin[j, m__max, 4] / 2,
                "i_minus_1": Dmin[j, m__max, 5] / 2,
                "C": Dmin[j, m__max, 6] / 2,
                "cts": int(i_Emax[j, 1026]),
                "di_over_2pi_percent": float(D_i_c_tot),
                "note": "Σ ∆i",
                "Polynom":"",
            })
            
            rows.append({
                "particle": "",
                "row_type": "total",
                "theory": "",
                "E": "",
                "E_m": Emin[j, 1028],
                "total": int(i_Emin[j, 0]),
                "i4": Dmin[j, m__min, 0] / 2,
                "i3": Dmin[j, m__min, 1] / 2,
                "i2": Dmin[j, m__min, 2] / 2,
                "i1": Dmin[j, m__min, 3] / 2,
                "i0": Dmin[j, m__min, 4] / 2,
                "i_minus_1": Dmin[j, m__min, 5] / 2,
                "C": Dmin[j, m__min, 6] / 2,
                "cts": int(i_Emin[j, 1026]),
                "di_over_2pi_percent": float(D_i_c_tot),
                "note": "Σ ∆i",
                "Polynom":"",
            })

    fieldnames = [
        "particle",
        "row_type",
        "theory",
        "E",
        "E_m",
        "total",
        "i4",
        "i3",
        "i2",
        "i1",
        "i0",
        "i_minus_1",
        "C",
        "cts",
        "di_over_2pi_percent",
        "note",
        "Polynom",
    ]

    if rows:
        path = os.path.join(RESULTS_DIR, f"{base_name}.csv")
        with open(path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)

    return labels
