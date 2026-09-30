# rare-disease-acmg-benchmark
#
# PYTHON must have numpy and matplotlib. A local virtualenv is the
# easiest route and is gitignored:
#     python3 -m venv .venv && .venv/bin/pip install numpy matplotlib
PYTHON ?= $(shell test -x .venv/bin/python && echo .venv/bin/python || echo python3)

.PHONY: help test figures analyse all

help:
	@echo "make test     - engine unit tests and results-integrity checks"
	@echo "make analyse  - Arm A statistics and Arm B stratifications"
	@echo "make figures  - regenerate results/figures/*.png"
	@echo "make all      - test, analyse, figures"
	@echo ""
	@echo "All of the above run on a fresh clone with no network and no external"
	@echo "resources. The numbered pipeline steps in scripts/ do not - see README."

test:
	$(PYTHON) tests/test_acmgx.py
	$(PYTHON) tests/test_results_integrity.py

analyse:
	$(PYTHON) scripts/14_analyse_arma.py
	$(PYTHON) scripts/17_armb_strata.py

figures:
	$(PYTHON) scripts/15_figures.py

all: test analyse figures
