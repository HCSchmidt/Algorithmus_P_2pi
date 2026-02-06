import numpy as np
import matplotlib.pyplot as plt
import csv
from .energy_model import EnergyModel
from .engine import ScanEngine


def run_sensitivity(
    *,
    sector,
    preset,
    particles,
    eps=1e-3,
    steps=21,
):
    """
    Sweep the base around 2π and measure sensitivity of the scan.
    """
    base_scales = np.linspace(1.0 - eps, 1.0 + eps, steps)
    rows = []

    for scale in base_scales:
        model = EnergyModel(base_scale=scale)
        engine = ScanEngine(
            preset=preset,
            model=model,
        )

        out = engine.run(particles)

        rows.append({
            "base_scale": scale,
            "possible_ET": out.possible_ET,
            "real_ET": out.real_ET,
            "hit_ratio": (
                out.real_ET / out.possible_ET
                if out.possible_ET > 0 else 0.0
            ),
        })

    return rows


def write_sensitivity_csv(path, rows):
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "base_scale",
                "possible_ET",
                "real_ET",
                "hit_ratio",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)



def plot_sensitivity(rows, out_png):
    x = [r["base_scale"] for r in rows]
    y_hits = [r["real_ET"] for r in rows]
    y_ratio = [r["hit_ratio"] for r in rows]

    fig, ax1 = plt.subplots(figsize=(8, 5))

    ax1.plot(x, y_hits, marker="o", label="real_ET")
    ax1.set_xlabel("Base scale (relative to 2π)")
    ax1.set_ylabel("real_ET")
    ax1.grid(True)

    ax2 = ax1.twinx()
    ax2.plot(x, y_ratio, linestyle="--", label="hit_ratio")
    ax2.set_ylabel("hit_ratio")

    fig.suptitle("Sensitivity to base variation")

    fig.tight_layout()
    fig.savefig(out_png, dpi=120)