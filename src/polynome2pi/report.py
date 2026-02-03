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

    with open(file_path, "w", encoding="utf8") as f:
        print("..............   plotting   .................")
        # (i_T / i_T1 prints are done in main because it has those values)

        print(
            "{0:10}{1:5}{2:20}{3:16}{4:8}{5:5}{6:5}{7:5}{8:5}{9:5}{10:5}{11:5}".format(
                "particle",
                "",
                "",
                "     theory: E",
                "   total ",
                "  i4",
                "  i3",
                "  i2",
                "  i1",
                "  i0",
                " i-1",
                "  C",
            ),
            file=f,
        )

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
                if flag == 0:
                    print(Obj[j][0], file=f)
                    flag += 1

                E_mean = (Emax[j, m] + Emin[j, m]) / 2
                Di_E = i_Emax[j, m] - i_Emin[j, m] + 1
                i_Emax[j, 0] += abs(Di_E)
                Cnt_ = Cnt[m]
                D_i_c = round((abs(Di_E)) * 100 / Cnt_, 5)
                i_Emax[j, 516] += Cnt_
                D_i_N[j] = float(i_Emax[j, 0]) * 100 / Cnt[m]

                Emax[j, m] = float(str(Emax[j, m])[:p_g])
                Emin[j, m] = float(str(Emin[j, m])[:p_g])
                E_mean = float(str(E_mean)[:p_g])

                print(
                    "{0:10}{1:5}{2:20}{3:16}{4:8}{5:5}{6:5}{7:5}{8:5}{9:5}{10:5}{11:5}".format(
                        "",
                        "max   ",
                        m_max,
                        Emax[j, m],
                        i_Emax[j, m],
                        Dmax[j, m, 0] / 2,
                        Dmax[j, m, 1] / 2,
                        Dmax[j, m, 2] / 2,
                        Dmax[j, m, 3] / 2,
                        Dmax[j, m, 4] / 2,
                        Dmax[j, m, 5] / 2,
                        Dmax[j, m, 6] / 2,
                    ),
                    file=f,
                )
                print(
                    "{0:10}{1:5}{2:20}{3:16}{4:8}{5:9}{6:7}{7:2}{8:12}{9:7}{10:2}".format(
                        "", "mean  ", Obj[j][1], E_mean, Di_E, "", "", "", "", "", ""
                    ),
                    file=f,
                )
                print(
                    "{0:10}{1:5}{2:20}{3:16}{4:8}{5:5}{6:5}{7:5}{8:5}{9:5}{10:5}{11:5}".format(
                        "",
                        "min   ",
                        m_min,
                        Emin[j, m],
                        i_Emin[j, m],
                        Dmin[j, m, 0] / 2,
                        Dmin[j, m, 1] / 2,
                        Dmin[j, m, 2] / 2,
                        Dmin[j, m, 3] / 2,
                        Dmin[j, m, 4] / 2,
                        Dmin[j, m, 5] / 2,
                        Dmin[j, m, 6] / 2,
                    ),
                    file=f,
                )
                print(
                    "{0:10}{1:5}{2:20}{3:15}{4:8}{5:5}{6:8}{7:9}{8:10}{9:2}".format(
                        "",
                        "",
                        "",
                        "         ∆ abs(i)",
                        abs(Di_E),
                        "  Cts",
                        Cnt_,
                        "  ∆i/(2pi)",
                        D_i_c,
                        " %",
                    ),
                    file=f,
                )

                # Restore labels in the plot
                if flag == 1:
                    labels.append((sector, j, i_Emin[j, m], E_mean, Obj[j][9], colors[j]))
                    flag = 2

            if j > 1 and m == 511 and i_Emax[j, 516] > 0:
                D_i_c = round(float(i_Emax[j, 0]) * 100 / i_Emax[j, 516], 5)
                D = str(D_i_c)
                D_i_c_[j] = "".rjust(10 - len(D), " ") + D + " %"
                print(
                    "{0:10}{1:5}{2:20}{3:15}{4:8}{5:5}{6:8}{7:9}{8:10}{9:2}".format(
                        "",
                        "total",
                        "",
                        "             Σ ∆i",
                        i_Emax[j, 0],
                        "  Cts",
                        i_Emax[j, 516],
                        "  ∆i/(2pi)",
                        D_i_c,
                        " %",
                    ),
                    file=f,
                )

            if flag == 0:
                print(
                    "{0:10}{1:5}{2:20}{3:16}{4:8}".format(
                        Obj[j][0], "mean  ", Obj[j][1], "  ", " only with i4 > 1"
                    ),
                    file=f,
                )
    return labels