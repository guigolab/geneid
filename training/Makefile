.PHONY: install test test-all lint fmt clean

install:
	python3 -m pip install -e ".[dev]"

# fast unit suite (skips the opt-in integration round-trip over real params)
test:
	python3 -m pytest -q -m "not integration"

# everything, including byte-identical round-trip over every param/*.param in the repo
test-all:
	python3 -m pytest -q

lint:
	python3 -m ruff check .

fmt:
	python3 -m ruff format .

clean:
	rm -rf build dist src/*.egg-info .pytest_cache
	find . -name __pycache__ -type d -prune -exec rm -rf {} +
