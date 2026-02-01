.PHONY: venv install run build test clean

venv:
	python3 -m venv .venv

install: venv
	. .venv/bin/activate && pip install -U pip && pip install -e .

run:
	. .venv/bin/activate && polynome2pi

build:
	. .venv/bin/activate && pip install -U build && python -m build

test:
	. .venv/bin/activate && pip install -U pytest && pytest -q

clean:
	rm -rf build dist .pytest_cache *.egg-info
	find . -name "__pycache__" -type d -exec rm -rf {} +