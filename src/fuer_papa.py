"""Entrypoint for Helmut :)"""

from polynome2pi.main import main

# =========================
# CONFIGURATION START
# =========================

ARGS = [
    "--sector", "nucleon",
    # "--no-open",
]

# Examples:
# ARGS = ["--sector", "minimal"]
# ARGS = ["--sector", "nucleon"]
# ARGS = ["--sector", "heavy"]
# ARGS = ["--sector", "light"]

# =======================
# USER CONFIGURATION END
# =======================

if __name__ == "__main__":
    main(["--sector", "minimal"])
    main(["--sector", "light"])
    main(["--sector", "broad"])
    main(["--sector", "nucleon"])
    main(["--sector", "heavy"])
