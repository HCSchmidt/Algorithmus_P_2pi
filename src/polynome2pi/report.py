import csv
from .cli import ScanSector
from .constants import get_colors


def write_report(
    file_path,
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
    colors = get_colors()
    labels = []

    rows = []

    for j in range(1, 29):
        m_min = float(Obj[j][2]) + float(Obj[j][3])
        m_max = float(Obj[j][2]) + float(Obj[j][4])
        p_g = len(Obj[j][2])
        m_min = float(str(m_min)[:p_g])
        m_max = float(str(m_max)[:p_g])

        flag = 0

        for m in range(1, 512):
            if i_Emax[j, m] == 0:
                continue

            E_mean = (Emax[j, m] + Emin[j, m]) / 2
            Di_E = i_Emax[j, m] - i_Emin[j, m] + 1
            i_Emax[j, 0] += abs(Di_E)

            Cnt_ = Cnt[m]
            D_i_c = round((abs(Di_E)) * 100 / Cnt_, 5)
            i_Emax[j, 516] += Cnt_
            D_i_N[j] = float(i_Emax[j, 0]) * 100 / Cnt[m]

            Emax_val = float(str(Emax[j, m])[:p_g])
            Emin_val = float(str(Emin[j, m])[:p_g])
            E_mean_val = float(str(E_mean)[:p_g])

            row = {
                "particle": Obj[j][0],
                "particle_symbol": Obj[j][9],
                "theory_E": Obj[j][2],
                "E_min": Emin_val,
                "E_mean": E_mean_val,
                "E_max": Emax_val,
                "delta_i": Di_E,
                "counts": Cnt_,
                "delta_i_percent": D_i_c,
                "i4": Dmax[j, m, 0] / 2,
                "i3": Dmax[j, m, 1] / 2,
                "i2": Dmax[j, m, 2] / 2,
                "i1": Dmax[j, m, 3] / 2,
                "i0": Dmax[j, m, 4] / 2,
                "i_minus_1": Dmax[j, m, 5] / 2,
                "C": Dmax[j, m, 6] / 2,
            }

            rows.append(row)

            # restore plot labels (unchanged behavior)
            if flag == 0:
                labels.append(
                    (sector, j, i_Emin[j, m], E_mean_val, Obj[j][9], colors[j])
                )
                flag = 1

        if j > 1 and i_Emax[j, 516] > 0:
            D_i_c_tot = round(float(i_Emax[j, 0]) * 100 / i_Emax[j, 516], 5)
            D_i_c_[j] = f"{D_i_c_tot} %"

    # ---- write CSV ----
    if rows:
        with open(file_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=rows[0].keys())
            writer.writeheader()
            writer.writerows(rows)

    return labels