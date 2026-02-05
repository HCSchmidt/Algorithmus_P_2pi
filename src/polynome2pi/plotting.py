from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Tuple
import webbrowser

import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

from .particles import Particle


colors =  [
        "#000000", "#F60000", "#05FB4F", "#CFCF00", "#07FCE4", "#F700D2",
        "#00F73E", "#7BB91F", "#A9BF06", "#789E20", "#CF00B7", "#EC61A9",
        "#FA9805", "#4200F6", "#495999", "#B91F50", "#CB4088", "#4D8E2F",
        "#499999", "#146108",
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
            plt.scatter(xs, ys, s=40, c=colors[idx % len(colors)], marker=".", linewidths=0, label=particles[key].name)

    if unmatched_segments:
        lc = LineCollection(unmatched_segments, colors="#C0BCBC", linewidths=3)
        plt.gca().add_collection(lc)

    # keep legend lightweight
    plt.legend(loc="best", fontsize=8, frameon=False)

    fig.savefig(out_png, dpi=120)
    plt.close(fig)


def open_file(path: Path) -> None:
    # OS-independent: use default app / browser
    webbrowser.open(path.resolve().as_uri())