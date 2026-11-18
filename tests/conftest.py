import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

# Makes the frozen reference implementation (tests/reference_*.py) importable.
sys.path.insert(0, str(Path(__file__).parent))
