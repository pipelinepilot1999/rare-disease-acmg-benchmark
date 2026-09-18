# rare-disease-acmg-benchmark
#
# PYTHON must have numpy, matplotlib and reportlab. A local virtualenv is the
# easiest route and is gitignored:
#     python3 -m venv .venv && .venv/bin/pip install numpy matplotlib reportlab
PYTHON ?= $(shell test -x .venv/bin/python && echo .venv/bin/python || echo python3)

.PHONY: help test guide figures analyse all clean-guide

help:
	@echo "make test     - engine unit tests and results-integrity checks"
	@echo "make analyse  - Arm A statistics and Arm B stratifications"
	@echo "make figures  - regenerate results/figures/*.png"
	@echo "make guide    - rebuild docs/interview-guide.pdf from results/"
	@echo "make all      - test, analyse, figures, guide"
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

guide: docs/interview-guide.pdf

docs/interview-guide.pdf: scripts/build_interview_guide.py \
                          results/ceiling.tsv results/arm_b_calls.tsv \
                          results/arm_a_ranks.tsv results/arm_a_summary.json \
                          results/criterion_metrics.tsv results/arm_b_strata.tsv \
                          src/acmgx/combine.py src/acmgx/assessability.py
	$(PYTHON) scripts/build_interview_guide.py

all: test analyse figures guide

clean-guide:
	rm -f docs/interview-guide.pdf
