# Makefile for polynome2pi
#
# IMPORTANT CONCEPTS:
#
# make install
#   - Creates a virtual environment (.venv)
#   - Installs this project INTO that environment (editable mode)
#   - This step creates the `polynome2pi` command
#   - Use this when you want to RUN the program
#
# make build
#   - Creates distributable files (wheel + source archive) in dist/
#   - Used for packaging or sharing, NOT for running locally
#   - Does NOT create or update the `polynome2pi` command
#
# In short:
#   install = for running
#   build   = for packaging

.PHONY: venv install run build test clean

# Create a local Python virtual environment in .venv/
venv:
	python3 -m venv .venv

# Install the project into the virtual environment (creates `polynome2pi` command)
install: venv
	. .venv/bin/activate && pip install -U pip && pip install -e .

# Run the installed program using the virtual environment
run:
	. .venv/bin/activate && polynome2pi

# Build distributable packages (wheel + sdist) into dist/
build:
	. .venv/bin/activate && pip install -U build && python -m build

# Run tests; installs pytest if needed)
test:
	. .venv/bin/activate && pip install -U pytest && pytest -q

# Remove build artifacts and caches
clean:
	rm -rf build dist .pytest_cache *.egg-info
	find . -name "__pycache__" -type d -exec rm -rf {} +