# Slides documenting the approximation error checks of the accuracy branch of
# Dynare.
#
#   make              the slides, from the results versioned in mod/results
#                     (LaTeX only, no MATLAB)
#   make results      run every example (MATLAB + Dynare, about forty minutes)
#   make <example>    run one example: growth linear burnside dhm ep_burnside
#                     irr bs
#   make dynare       one-time setup of a worktree of the accuracy branch
#   make clean        LaTeX and Dynare temporary files
#   make distclean    also the PDF, the stamps and the .mat files of the results

DYNARE_REPO      ?= $(HOME)/works/Dynare/dynare
DYNARE_ROOT      ?= $(HOME)/works/Dynare/accuracy
DYNARE_BUILD_DIR ?= $(DYNARE_ROOT)/build-matlab
MATLAB_ROOT      ?= $(HOME)/matlab/R2025b
MATLAB           ?= $(MATLAB_ROOT)/bin/matlab

export DYNARE_ROOT DYNARE_BUILD_DIR

EXAMPLES = growth linear burnside dhm ep_burnside irr bs
STAMPS   = $(EXAMPLES:%=mod/results/%.stamp)

RUN = cd mod && $(MATLAB) -batch "run_example('$(1)')"

.PHONY: all results slides dynare clean distclean $(EXAMPLES)

all: slides

slides: tex/accuracy.pdf

# The slides only read the versioned results: building them never runs an
# example, whatever the timestamps left by a checkout.
RESULTS = $(wildcard mod/results/*.tex mod/results/*.csv) $(wildcard mod/results/dynare.hash)

tex/accuracy.pdf: tex/accuracy.tex $(RESULTS)
	cd tex && rubber --pdf --unsafe accuracy.tex

# Running the examples. A stamp records that an example was run with its current
# sources; the stamps are not versioned, so after a checkout each example runs
# once.
results: $(STAMPS)

$(EXAMPLES): %: mod/results/%.stamp

mod/results/%.stamp: mod/%.mod mod/run_example.m $(wildcard mod/export/*.m)
	$(call RUN,$*)

mod/results/burnside.stamp mod/results/dhm.stamp: mod/burnside_exact_policy.m

IRR_SOURCES = mod/occbin_irr.mod mod/ep_irr.mod mod/sep_irr.mod mod/sep_irr_h2.mod mod/sep_irr_h4.mod \
              mod/compare_irr.mod mod/sep_irr.inc mod/sweep_irr.mod mod/recalibrate.m \
              mod/irr_declarations.inc mod/irr_calibration.inc mod/irr_min_model.inc

mod/results/irr.stamp: $(IRR_SOURCES) mod/run_example.m $(wildcard mod/export/*.m)
	$(call RUN,irr)

BS_SOURCES = mod/occbin_bs.mod mod/ep_bs.mod mod/sep_bs.mod mod/sep_bs_h2.mod mod/sep_bs_h4.mod \
             mod/sep_bs_o2.mod mod/sep_bs_o3.mod mod/compare_bs.mod mod/sep_bs.inc mod/bs_declarations.inc \
             mod/bs_calibration.inc mod/bs_min_model.inc mod/buffer_egm.m mod/buffer_policy.m \
             mod/buffer_exact_policy.m mod/score_rules.m

mod/results/bs.stamp: $(BS_SOURCES) mod/run_example.m $(wildcard mod/export/*.m)
	$(call RUN,bs)

# One-time setup: a worktree of the accuracy branch, built for MATLAB. The main
# checkout of Dynare is left alone.
dynare: $(DYNARE_BUILD_DIR)/preprocessor/src/dynare-preprocessor

$(DYNARE_ROOT)/meson.build:
	git -C $(DYNARE_REPO) worktree add $(DYNARE_ROOT) accuracy
	cd $(DYNARE_ROOT) && git submodule update --init --recursive

$(DYNARE_BUILD_DIR)/preprocessor/src/dynare-preprocessor: $(DYNARE_ROOT)/meson.build
	cd $(DYNARE_ROOT) && meson setup -Dbuild_for=matlab -Dmatlab_path=$(MATLAB_ROOT) \
	  --buildtype=debugoptimized $(DYNARE_BUILD_DIR)
	ninja -C $(DYNARE_BUILD_DIR)

# Dynare writes a +<name> package and a <name> folder for each .mod file run in
# mod/; every folder there is generated, but for export/ and results/.
clean:
	cd tex && rubber --clean --pdf accuracy.tex || true
	rm -f tex/git.hash tex/*.synctex.gz
	find mod -mindepth 1 -maxdepth 1 -type d ! -name export ! -name results -exec rm -rf {} +
	rm -f mod/*.log mod/ctx_*.mat

distclean: clean
	rm -f tex/accuracy.pdf mod/results/*.stamp mod/results/*.mat
