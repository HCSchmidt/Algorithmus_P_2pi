from __future__ import annotations

import webbrowser
from pathlib import Path


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


def open_image(path: Path) -> None:
    # OS-independent: use default app / browser
    webbrowser.open(path.resolve().as_uri())
