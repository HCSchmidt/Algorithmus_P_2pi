from __future__ import annotations
import numpy as np
from typing import Optional
from pathlib import Path
from typing import Dict, List, Tuple
from polynome2pi.engine import ScanOutputs
from polynome2pi.presets import ScanPreset  
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401  (needed for 3D)


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
    sector_name: str,
) -> None:
    out_png.parent.mkdir(parents=True, exist_ok=True)

    fig = plt.gcf()
    fig.set_size_inches(10, 6)

    title = (f"P(2π) scan: sector {sector_name}",)

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

def plot_match_grid_3d_scatter(
    outputs: ScanOutputs,
    preset: ScanPreset,
    particles: dict[str, Particle],
    path: Path,
    particle_key: Optional[str] = None,
    min_hits: int = 1,
) -> None:
    """
    3D scatter over the true grid (i4, i3, i2).

    Axes:
      x = i3
      y = i4
      z = i2

    Color encodes hit count per grid cell.
    Only cells with hit_count >= min_hits are drawn.

    If particle_key is None -> ANY particle grid
    else -> grid for that particle_key.
    """
    if particle_key is None:
        grid3 = outputs.any_match_grid
        title = f"3D Match-Grid (ANY): sector{preset.sector.value}"
    else:
        if particle_key not in outputs.particle_match_grids:
            raise KeyError(f"Unknown particle_key '{particle_key}'.")
        grid3 = outputs.particle_match_grids[particle_key]
        p = particles[particle_key]
        title = f"3D Match-Grid ({p.name}) – sector={preset.sector.value}"

    # find nonzero (or >= min_hits) cells
    mask = grid3 >= int(min_hits)
    if not np.any(mask):
        # still write an empty plot to make it obvious
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
        ax.set_title(title + f" (no cells >= {min_hits})")
        ax.set_xlabel("i3")
        ax.set_ylabel("i4")
        ax.set_zlabel("i2")
        fig.set_size_inches(10, 7)
        path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        return

    gi4_idx, gi3_idx, gi2_idx = np.where(mask)
    counts = grid3[gi4_idx, gi3_idx, gi2_idx].astype(float)

    # map grid indices -> actual coefficient values
    i4_vals = np.array(outputs.grid_i4_values)
    i3_vals = np.array(outputs.grid_i3_values)
    i2_vals = np.array(outputs.grid_i2_values)

    y_i4 = i4_vals[gi4_idx]
    x_i3 = i3_vals[gi3_idx]
    z_i2 = i2_vals[gi2_idx]

    # size scaling: keep readable across sectors
    # (sqrt to compress large outliers)
    sizes = 10.0 + 40.0 * np.sqrt(counts / counts.max())

    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    
    _draw_const_planes(ax, i4_vals, i3_vals, i2_vals, alpha=0.1, lw=0.5)

    sc = ax.scatter(
        x_i3,
        y_i4,
        z_i2,
        c=counts,
        s=sizes,
        marker="o",
        linewidths=0,
        alpha=0.9,
        cmap="viridis",
    )

    ax.set_title(title)
    ax.set_xlabel("i3")
    ax.set_ylabel("i4")
    ax.set_zlabel("i2")

    # optional: make ticks less crazy (only show actual scanned values)
    ax.set_xlim(min(outputs.grid_i3_values), max(outputs.grid_i3_values))
    ax.set_ylim(min(outputs.grid_i4_values), max(outputs.grid_i4_values))
    ax.set_zlim(min(outputs.grid_i2_values), max(outputs.grid_i2_values))


    cbar = fig.colorbar(sc, ax=ax, shrink=0.75, pad=0.1)
    cbar.set_label("hit count per (i4,i3,i2) cell")

    fig.set_size_inches(10, 7)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)

def _draw_const_planes(
    ax,
    i4_vals,
    i3_vals,
    i2_vals,
    alpha=0.15,
    lw=0.5,
):
    # i2 = const → Linien im (i3,i4)-Raster
    for z in i2_vals:
        for y in i4_vals:
            ax.plot(
                i3_vals,
                [y] * len(i3_vals),
                [z] * len(i3_vals),
                color="grey",
                alpha=alpha,
                linewidth=lw,
            )

    # i3 = const → Linien im (i2,i4)-Raster
    for x in i3_vals:
        for y in i4_vals:
            ax.plot(
                [x] * len(i2_vals),
                [y] * len(i2_vals),
                i2_vals,
                color="grey",
                alpha=alpha,
                linewidth=lw,
            )

    # i4 = const → Linien im (i2,i3)-Raster
    for y in i4_vals:
        for x in i3_vals:
            ax.plot(
                [x] * len(i2_vals),
                [y] * len(i2_vals),
                i2_vals,
                color="grey",
                alpha=alpha,
                linewidth=lw,
            )