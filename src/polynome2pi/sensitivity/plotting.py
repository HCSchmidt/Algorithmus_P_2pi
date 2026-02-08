from polynome2pi.sensitivity.SensitivityPoint import SensitivityPoint


import matplotlib.pyplot as plt


from pathlib import Path
from typing import List


def plot_sensitivity(
    points: List[SensitivityPoint],
    png_path: Path,
    *,
    title: str,
    y_left: str = "real_ET",
    y_right: str = "hit_ratio",
) -> None:
    """
    y_left:
      - "real_ET" (global)
      - "particle_hits" (only if particle_key set)
    """
    xs = [p.base_scale for p in points]

    if y_left == "real_ET":
        y1 = [p.real_ET for p in points]
        y1_label = "real_ET"
    elif y_left == "particle_hits":
        y1 = [p.particle_hits or 0 for p in points]
        y1_label = "particle_hits"
    else:
        raise ValueError(f"Unknown y_left: {y_left}")

    if y_right == "hit_ratio":
        y2 = [p.hit_ratio for p in points]
        y2_label = "hit_ratio"
    elif y_right == "particle_hit_ratio":
        y2 = [p.particle_hit_ratio or 0.0 for p in points]
        y2_label = "particle_hit_ratio"
    else:
        raise ValueError(f"Unknown y_right: {y_right}")

    png_path.parent.mkdir(parents=True, exist_ok=True)

    fig, ax1 = plt.subplots(figsize=(10, 6))
    ax1.plot(xs, y1, marker="o")
    ax1.set_xlabel("Base scale (relative to 2π)")
    ax1.set_ylabel(y1_label)
    ax1.grid(True)

    ax2 = ax1.twinx()
    ax2.plot(xs, y2, linestyle="--")
    ax2.set_ylabel(y2_label)

    plt.title(title)
    plt.tight_layout()
    plt.savefig(png_path, dpi=120)
    plt.close(fig)
