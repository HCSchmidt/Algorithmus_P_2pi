from __future__ import annotations

import webbrowser
from pathlib import Path


def open_image(path: Path) -> None:
    # OS-independent: use default app / browser
    webbrowser.open(path.resolve().as_uri())
