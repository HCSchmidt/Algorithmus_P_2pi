import os
import sys


def open_image(path: str) -> None:
    """Open an image file with the OS default viewer (cross-platform)."""
    if os.name == "nt":  # Windows
        os.startfile(path)
    elif sys.platform == "darwin":  # macOS
        os.system(f'open "{path}"')
    else:  # Linux, WSL, others
        os.system(f'xdg-open "{path}"')