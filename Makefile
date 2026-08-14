PYTHON ?= python

.PHONY: help figures figures-published supplementary all test verify notebook notebook-build clean

help:
	@echo "Every target below has an equivalent: python run_all.py <target>"
	@echo "(use that on Windows, where make is usually not installed)"
	@echo ""
	@echo "make figures            main figures 1-5 into figures/_out/"
	@echo "make supplementary      supplementary figures S1-S6"
	@echo "make all                both of the above"
	@echo "make figures-published  regenerate the submitted Fig 4B/4C (wrong control clock)"
	@echo "make test               run the test suite"
	@echo "make verify             run the test suite and print the published statistics"
	@echo "make notebook           launch the panel-by-panel walkthrough"
	@echo "make notebook-build     rebuild that notebook from tools/ and execute it"
	@echo "make clean              remove caches (figures/_out/ is committed, so it stays)"

figures:
	$(PYTHON) figures/figure1.py
	$(PYTHON) figures/figure2.py
	$(PYTHON) figures/figure3.py
	$(PYTHON) figures/figure4.py
	$(PYTHON) figures/figure5.py

supplementary:
	$(PYTHON) figures/supplementary/figS1_barcode_linear.py
	$(PYTHON) figures/supplementary/figS2_cluster_selection.py
	$(PYTHON) figures/supplementary/figS3_coclustering_support.py
	$(PYTHON) figures/supplementary/figS4_window_sweep.py
	$(PYTHON) figures/supplementary/figS5_eigenvalue_spectra.py
	$(PYTHON) figures/supplementary/figS6_global_coclustering.py

all: figures supplementary

# The submitted Figure 4B/4C, which puts control indices through the colonised
# clock. Kept so the difference can be seen rather than only described.
figures-published:
	$(PYTHON) figures/figure4.py --clock published

test:
	$(PYTHON) -m pytest tests/ -q

verify:
	$(PYTHON) -m pytest tests/ -v

notebook:
	jupyter lab notebooks/reproduce_figures.ipynb

# The notebook is generated, not hand-edited. Edit tools/build_notebook.py.
notebook-build:
	$(PYTHON) tools/build_notebook.py --execute

# figures/_out/ is committed, so it is not cleaned: `git checkout figures/_out`
# is the way back if a run leaves it dirty.
clean:
	rm -rf .pytest_cache
	find . -type d -name __pycache__ -prune -exec rm -rf {} +
	find . -type d -name .ipynb_checkpoints -prune -exec rm -rf {} +
