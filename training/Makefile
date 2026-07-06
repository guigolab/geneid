# geneid-train dev tooling. Requires Python >= 3.10 (see pyproject.toml).
#
# Quick start (from this directory):
#     make venv                    # if `python3` is already >= 3.10
#     make venv PYTHON=python3.12  # else point at a 3.10+ interpreter
#     make test
#
# `make venv` builds an isolated ./.venv; the other targets auto-use it, so
# there is nothing to activate. macOS ships 3.9 as `python3` -- if so, install a
# newer Python (e.g. `brew install python@3.12`) and pass PYTHON=python3.12.

PYTHON ?= python3
VENV   := .venv

# Prefer the project venv if it exists, else fall back to the active interpreter.
PY := $(shell [ -x $(VENV)/bin/python ] && echo $(VENV)/bin/python || echo $(PYTHON))

.PHONY: venv install test test-all lint fmt clean

# One-command setup: isolated venv + editable install with dev tools (pytest, ruff).
# No activation needed afterwards -- `make test`/`lint`/... pick up ./.venv.
venv:
	@$(PYTHON) -c 'import sys; raise SystemExit(0 if sys.version_info >= (3,10) else 1)' || { \
	  echo "ERROR: '$(PYTHON)' is $$($(PYTHON) -V 2>&1); geneid-train needs Python >= 3.10."; \
	  echo "       Install a newer Python (e.g. 'brew install python@3.12') and re-run:"; \
	  echo "           make venv PYTHON=python3.12"; \
	  exit 1; }
	$(PYTHON) -m venv $(VENV)
	$(VENV)/bin/python -m pip install --upgrade pip
	$(VENV)/bin/python -m pip install -e ".[dev]"
	@echo "Ready. Run 'make test', or './.venv/bin/geneid-train --help'."

# Editable install + dev deps into the CURRENT environment (use inside your own venv).
install:
	$(PY) -m pip install --upgrade pip
	$(PY) -m pip install -e ".[dev]"

# fast unit suite (skips the opt-in integration round-trip over real params)
test:
	$(PY) -m pytest -q -m "not integration"

# everything, including byte-identical round-trip over every param/*.param in the repo
test-all:
	$(PY) -m pytest -q

lint:
	$(PY) -m ruff check .

fmt:
	$(PY) -m ruff format .

clean:
	rm -rf build dist src/*.egg-info .pytest_cache
	find . -name __pycache__ -type d -prune -exec rm -rf {} +
