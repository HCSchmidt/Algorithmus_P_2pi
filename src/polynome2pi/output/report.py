from __future__ import annotations

import csv
from typing import Dict

from ..cli import ScanSector
from ..particles import Particle


CSV_HEADER = [
    "particle",
    "row_type",
    "theory",
    "E",
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
]


def _truncate_like_legacy(value: float, p_g: int) -> float:
    # legacy did: float(str(value)[:p_g])
    # (works for your examples: 10.0, 5.14, 0.995, etc.)
    return float(str(float(value))[:p_g])


def write_report_csv(
    file_path,
    sector: ScanSector,
    particles: Dict[str, Particle],
    results: Dict[str, object],  # ParticleData-like
    Cnt,
):
    """
    Writes LEGACY CSV format (the one you showed under 'before').

    Assumes results[p.key] has:
      - Emax[m], Emin[m], i_Emax[m], i_Emin[m]
      - Dmax[m], Dmin[m] each with 7 entries: i4,i3,i2,i1,i0,i-1,C (raw ints)
      - m_values() -> range(1,512)
    """

    with open(file_path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(CSV_HEADER)

        for key, p in particles.items():
            pdata = results.get(key)

            # particle header row (legacy)
            w.writerow([p.name, "particle"] + [""] * (len(CSV_HEADER) - 2))

            # If no pdata or no matches -> legacy "only with i4 > 1"
            if pdata is None:
                w.writerow([p.name, "mean", p.theory_str] + [""] * (len(CSV_HEADER) - 4) + ["only with i4 > 1"])
                continue

            any_block = False
            total_delta_abs = 0
            total_cts = 0

            # legacy precision control
            p_g = len(str(p.theory_E))

            for m in pdata.m_values():
                if pdata.i_Emax[m] == 0:
                    continue

                any_block = True

                emax = float(pdata.Emax[m])
                emin = float(pdata.Emin[m]) if float(pdata.Emin[m]) != 0 else emax
                eme  = (emax + emin) / 2

                # legacy min/max theory bounds
                m_min = _truncate_like_legacy(p.theory_E + p.sd_minus, p_g)
                m_max = _truncate_like_legacy(p.theory_E + p.sd_plus, p_g)

                # legacy delta i
                di = int(pdata.i_Emax[m] - pdata.i_Emin[m] + 1)
                di_abs = abs(di)

                # counts for this m (legacy uses Cnt[m])
                cts = int(Cnt[m]) if int(Cnt[m]) != 0 else 0
                di_pct = round(di_abs * 100 / cts, 5) if cts else ""

                total_delta_abs += di_abs
                total_cts += cts

                # unpack D vectors (raw ints), legacy divides by 2 in output
                dmax = pdata.Dmax[m]
                dmin = pdata.Dmin[m]

                def _fmt_half(x):
                    return float(x) / 2

                # max row: theory=m_max, E=Emax, total=i_Emax
                w.writerow([
                    "", "max", m_max, _truncate_like_legacy(emax, p_g), int(pdata.i_Emax[m]),
                    _fmt_half(dmax[0]), _fmt_half(dmax[1]), _fmt_half(dmax[2]),
                    _fmt_half(dmax[3]), _fmt_half(dmax[4]), _fmt_half(dmax[5]),
                    _fmt_half(dmax[6]),  # C/2
                    "", "", ""
                ])

                # mean row: theory=theory_str, E=E_mean, total=Di_E (can be negative)
                w.writerow([
                    "", "mean", p.theory_str, _truncate_like_legacy(eme, p_g), di,
                    "", "", "", "", "", "", "",
                    "", "", ""
                ])

                # min row: theory=m_min, E=Emin, total=i_Emin
                w.writerow([
                    "", "min", m_min, _truncate_like_legacy(emin, p_g), int(pdata.i_Emin[m]),
                    _fmt_half(dmin[0]), _fmt_half(dmin[1]), _fmt_half(dmin[2]),
                    _fmt_half(dmin[3]), _fmt_half(dmin[4]), _fmt_half(dmin[5]),
                    _fmt_half(dmin[6]),  # C/2
                    "", "", ""
                ])

                # delta row: total=abs(Di_E), cts, di_over_2pi_percent, note="∆ abs(i)"
                w.writerow([
                    "", "delta", "", "", di_abs,
                    "", "", "", "", "", "", "",
                    cts, di_pct, "∆ abs(i)"
                ])

            # if nothing matched, legacy prints "only with i4 > 1"
            if not any_block:
                w.writerow([p.name, "mean", p.theory_str] + [""] * (len(CSV_HEADER) - 4) + ["only with i4 > 1"])
                continue

            # total row for that particle (legacy note Σ ∆i)
            total_pct = round(total_delta_abs * 100 / total_cts, 5) if total_cts else ""
            w.writerow([
                "", "total", "", "", total_delta_abs,
                "", "", "", "", "", "", "",
                total_cts, total_pct, "Σ ∆i"
            ])