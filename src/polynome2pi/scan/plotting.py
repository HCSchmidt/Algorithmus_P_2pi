from __future__ import annotations

from typing import Optional
from pathlib import Path
from typing import Dict, List, Tuple
from polynome2pi.engine import ScanOutputs
from polynome2pi.presets import ScanPreset  
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

from ..particles import Particle


colors = [
    "#000000",
    "#F60000",
    "#05FB4F",
    "#CFCF00",
    "#07FCE4",
    "#F700D2",
    "#00F73E",
    "#7BB91F",
    "#A9BF06",
    "#789E20",
    "#CF00B7",
    "#EC61A9",
    "#FA9805",
    "#4200F6",
    "#495999",
    "#B91F50",
    "#CB4088",
    "#4D8E2F",
    "#499999",
    "#146108",
]


def plot_scan(
    out_png: Path,
    particles: Dict[str, Particle],
    matched_points: Dict[str, Tuple[List[int], List[float]]],
    unmatched_segments: List[Tuple[Tuple[int, float], Tuple[int, float]]],
    title: str,
) -> None:
    out_png.parent.mkdir(parents=True, exist_ok=True)

    fig = plt.gcf()
    fig.set_size_inches(10, 6)

    plt.title(title)
    plt.xlabel("Scan index N")
    plt.ylabel("Energy in $m_e$")

    keys = list(particles.keys())

    for idx, key in enumerate(keys):
        xs, ys = matched_points.get(key, ([], []))
        if xs:
            plt.scatter(
                xs,
                ys,
                s=40,
                c=colors[idx % len(colors)],
                marker=".",
                linewidths=0,
                label=particles[key].name,
            )

    if unmatched_segments:
        lc = LineCollection(unmatched_segments, colors="#C0BCBC", linewidths=3)
        plt.gca().add_collection(lc)

    # keep legend lightweight
    plt.legend(loc="best", fontsize=8, frameon=False)

    fig.savefig(out_png, dpi=120)
    plt.close(fig)


def plot_match_grid(
    outputs: ScanOutputs,
    preset: ScanPreset,
    particles: dict[str, Particle],
    path: Path,
    particle_key: Optional[str] = None,
) -> None:
    """
    Heatmap grid over (i4, i3), summing over i2:
      - if particle_key is None: counts where ANY particle matched
      - else: counts where that particle matched

    x-axis: i3
    y-axis: i4
    value: sum over i2 of hit-counts
    """
    if particle_key is None:
        grid3 = outputs.any_match_grid
        title = f"Match-Grid (ANY) – sector={preset.sector.value}"
    else:
        if particle_key not in outputs.particle_match_grids:
            raise KeyError(f"Unknown particle_key '{particle_key}'.")
        grid3 = outputs.particle_match_grids[particle_key]
        p = particles[particle_key]
        title = f"Match-Grid ({p.name}) – sector={preset.sector.value}"

    # reduce 3D -> 2D by summing over i2
    grid2 = grid3.sum(axis=2)  # shape (n_i4, n_i3)

    i4_vals = outputs.grid_i4_values
    i3_vals = outputs.grid_i3_values

    fig = plt.figure()
    ax = plt.gca()

    im = ax.imshow(
        grid2,
        origin="lower",
        aspect="auto",
        interpolation="nearest",
    )

    ax.set_title(title)
    ax.set_xlabel("i3")
    ax.set_ylabel("i4")

    ax.set_xticks(range(len(i3_vals)))
    ax.set_xticklabels(i3_vals)

    ax.set_yticks(range(len(i4_vals)))
    ax.set_yticklabels(i4_vals)

    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label("hit count (summed over i2)")

    fig.set_size_inches(10, 6)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)