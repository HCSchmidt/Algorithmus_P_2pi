"""Entrypoint for Helmut :)"""

from polynome2pi.main import main

# =========================
# CONFIGURATION START
# =========================

ARGS = [
    "--sector", "nucleon",
    "--no-show",
]

# Examples:
# ARGS = ["--sector", "minimal"]
# ARGS = ["--sector", "heavy", "--base-name", "tau_scan"]
# ARGS = ["--sector", "light"]

# =======================
# USER CONFIGURATION END
# =======================

if __name__ == "__main__":
    print("Running with arguments:", ARGS)
    main(ARGS)
