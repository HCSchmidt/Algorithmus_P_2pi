"""
Run every evaluation and write CSVs + figures to results/evaluation/.

    python -m polynome2pi.evaluation.run_all [--bases 200]
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from polynome2pi.evaluation.nulltest import (
    base_null,
    correction_series,
    evaluate_targets,
    holdout_targets,
    pdg_electron,
    random_targets,
    spectrum,
    table_targets,
    tuning_capacity,
)
from polynome2pi.evaluation.paper_model import coefficient_grid

BLUE, ORANGE, MUTED, INK = "#2a78d6", "#eb6834", "#c3c2b7", "#0b0b0b"
SHARP_LAMBDA = 0.1  # windows narrower than ~1/10 of the level spacing


def _style(ax):
    ax.grid(True, color="#e6e5e0", linewidth=0.6)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def run(out: Path, n_bases: int, n_random: int) -> None:
    out.mkdir(parents=True, exist_ok=True)
    grid = coefficient_grid()
    spec = spectrum(grid=grid)
    targets = table_targets() + [pdg_electron()] + holdout_targets()

    # 1. Every particle under 2π: hits vs chance expectation
    df = pd.DataFrame(evaluate_targets(spec, targets))
    df.to_csv(out / "particles_2pi.csv", index=False)
    print(df[["key", "group", "hits", "lam", "best_rel", "in_range"]].to_string())

    fig, ax = plt.subplots(figsize=(7, 5))
    d = df[df.in_range & (df.hits > 0)]
    lim = [1e-5, 1e5]
    ax.plot(lim, lim, color=MUTED, linewidth=1.5, label="chance: hits = λ")
    for grp, color in (("table", BLUE), ("holdout", ORANGE)):
        g = d[d.group == grp]
        ax.scatter(g.lam, g.hits, s=40, color=color, edgecolor="white", linewidth=1.5, zorder=3,
                   label={"table": "paper's particles", "holdout": "held-out particles"}[grp])
    for i, (_, r) in enumerate(d[d.lam < SHARP_LAMBDA].sort_values("lam").iterrows()):
        ax.annotate(r.key, (r.lam, r.hits), textcoords="offset points", xytext=(0, 10 if i % 2 == 0 else -16),
                    ha="center", fontsize=9, color=INK)
    ax.set(xscale="log", yscale="log", xlim=lim, ylim=[0.5, 1e5],
           xlabel="λ = expected chance hits (density × window width)", ylabel="observed hits in window",
           title="Hits track chance except for three fitted particles")
    _style(ax)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(out / "hits_vs_chance.png", dpi=150)
    plt.close(fig)

    # 2. Random targets with the same relative widths as the table
    tbl = df[(df.group == "table") & df.in_range & (df.width > 0)]
    rel_w = (tbl.width / tbl.mass).to_numpy()
    rnd = pd.DataFrame(evaluate_targets(spec, random_targets(spec, rel_w, n_random, 1.0, 3700.0)))
    summary = [
        f"random targets: n={len(rnd)} hit_rate={rnd.hit.mean():.3f} expected={rnd.p_chance.mean():.3f}",
        f"table particles in range: n={len(tbl)} hit_rate={tbl.hit.mean():.3f} expected={tbl.p_chance.mean():.3f}",
    ]

    # 3. Alternative bases (π := base/2 everywhere), sharp particles only
    sharp_keys = set(df[(df.lam < SHARP_LAMBDA) & df.in_range].key)
    sharp = [t for t in targets if t.key in sharp_keys]
    rng = np.random.default_rng(1)
    named = {"6": 6.0, "2e": 2 * math.e, "π+3": math.pi + 3, "√40": math.sqrt(40), "6.5": 6.5, "7": 7.0}
    bases = np.r_[2 * math.pi, list(named.values()), rng.uniform(5.5, 7.0, n_bases)]
    bn = pd.DataFrame(base_null(sharp, bases, grid))
    bn.to_csv(out / "base_null.csv", index=False)
    per_base = bn.groupby("base").agg(sharp_hits=("hit", "sum"), med_rel=("best_rel", "median"))
    tpi = per_base.loc[2 * math.pi]
    others = per_base.drop(index=2 * math.pi)
    summary += [
        f"base null: sharp targets={sorted(sharp_keys)}",
        f"  2π: sharp_hits={int(tpi.sharp_hits)} median best_rel={tpi.med_rel:.2e}",
        f"  other bases (n={len(others)}): mean sharp_hits={others.sharp_hits.mean():.3f}, "
        f"max={int(others.sharp_hits.max())}, p(hits>=2π)={(others.sharp_hits >= tpi.sharp_hits).mean():.4f}",
    ]
    for name, b in named.items():
        summary.append(f"  base {name}={b:.5f}: sharp_hits={int(per_base.loc[b].sharp_hits)}")

    fig, ax = plt.subplots(figsize=(7, 4.5))
    for key, color in (("proton", BLUE), ("neutron", ORANGE)):
        sub = bn[(bn.key == key) & (bn.base != 2 * math.pi)]
        ax.hist(np.log10(sub.best_rel), bins=np.arange(-11.5, -5, 0.2), histtype="step", linewidth=2,
                color=color, label=f"{key}: {len(sub)} other bases")
        two = bn[(bn.key == key) & (bn.base == 2 * math.pi)].best_rel.iloc[0]
        ax.axvline(math.log10(two), color=color, linewidth=2, linestyle="--", label=f"{key}: 2π")
    ax.set(xlabel="log10 relative distance of closest grid energy", ylabel="number of bases",
           title="Only 2π lands on proton and neutron (the constants were tuned there)")
    _style(ax)
    ax.legend(frameon=False, loc="upper left")
    fig.tight_layout()
    fig.savefig(out / "base_null.png", dpi=150)
    plt.close(fig)

    # 4. Tuning capacity: grid energy + one Σ a_k π^-k constant (k=7..14, a_k ∈ {-2..2})
    corr = correction_series()
    proton = next(t for t in targets if t.key == "proton")
    width = proton.hi - proton.lo
    cap = []
    for b_label, b in (("2π", 2 * math.pi), ("6.5", 6.5), ("random", float(rng.uniform(5.5, 7.0)))):
        spec_b = spec if b == 2 * math.pi else spectrum(b, grid)
        corr_b = corr if b == 2 * math.pi else correction_series(pi=b / 2)
        masses = rng.uniform(1800, 1900, 200)
        fake = [type(proton)(f"f{i}", m, m - width / 2, m + width / 2, "random") for i, m in enumerate(masses)]
        ok = np.mean([tuning_capacity(spec_b, corr_b, t) for t in fake])
        cap.append({"base": b_label, "value": b, "fraction_hit_at_proton_precision": ok})
        summary.append(f"tuning capacity base={b_label}: {ok:.3f} of 200 random targets hit at proton precision")
    pd.DataFrame(cap).to_csv(out / "tuning_capacity.csv", index=False)

    (out / "summary.txt").write_text("\n".join(summary) + "\n")
    print("\n".join(summary))


def main(argv=None) -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--bases", type=int, default=200)
    p.add_argument("--random-targets", type=int, default=20000)
    p.add_argument("--out", type=Path, default=Path("results/evaluation"))
    a = p.parse_args(argv)
    run(a.out, a.bases, a.random_targets)


if __name__ == "__main__":
    main()
