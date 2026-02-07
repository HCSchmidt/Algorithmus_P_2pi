from __future__ import annotations

import os
import sys
import subprocess
from pathlib import Path


# ---------------------------
# Cross-platform image opening
# ---------------------------

def open_image(path: str | Path) -> None:
    """
    Open an image file with the default system viewer (cross-platform).

    Works on:
      - macOS
      - Linux
      - Windows
    """
    path = str(path)

    try:
        if sys.platform.startswith("darwin"):
            subprocess.run(["open", path], check=False)
        elif os.name == "nt":
            os.startfile(path)  # type: ignore[attr-defined]
        elif os.name == "posix":
            subprocess.run(["xdg-open", path], check=False)
    except Exception:
        # Fail silently: opening the image is convenience only
        pass


# ---------------------------
# Small numeric helpers
# ---------------------------

def safe_div(numerator: float, denominator: float, default: float = 0.0) -> float:
    """
    Safe division helper to avoid ZeroDivisionError.
    """
    return numerator / denominator if denominator else default


def truncate_like_legacy(value: float, digits: int) -> float:
    """
    Truncate (not round!) a float to a given number of significant digits,
    mimicking the legacy string-slice behaviour.
    """
    s = f"{value:.20g}"
    return float(s[:digits]) if len(s) > digits else float(s)