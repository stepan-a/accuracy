# Accuracy of solutions in Dynare — slides

[![slides](https://github.com/stepan-a/accuracy/actions/workflows/slides.yml/badge.svg?branch=master)](https://github.com/stepan-a/accuracy/actions/workflows/slides.yml)

Beamer slides documenting the approximation error checks of the `accuracy`
branch of Dynare (`matlab/+accuracy/`): residuals and their normalizations,
quadrature, evaluation points, the Den Haan and Marcet test, and comparisons of
perturbation, extended path and OccBin solutions, illustrated by examples that
run.

**[Download the slides (PDF)](https://github.com/stepan-a/accuracy/releases/latest/download/accuracy.pdf)**,
built by the continuous integration from the last commit on `master`.

## Layout

- `tex/accuracy.tex` — the slides (English).
- `mod/` — the examples, adapted from `tests/accuracy/` of the branch, and their
  helpers. Each example writes the numbers quoted in the slides to
  `mod/results/` with the writers of `mod/export/` (LaTeX macros, booktabs
  tables, CSV files for pgfplots).
- `mod/results/` — the results read by the slides. They are versioned, so that
  the slides build without running the examples; `dynare.hash` records the
  commit of the branch that produced them, shown in the footer of the slides.
- `Makefile` — builds the slides and runs the examples.
- `.github/workflows/slides.yml` — builds the slides and publishes the PDF.

## Building the slides

```sh
make            # tex/accuracy.pdf, from the versioned results
```

This needs LaTeX only (`rubber`, `beamer`, `pgfplots`, `ccicons`), not MATLAB.
The slides call `git` through `\write18` to print the commit in their footer,
hence `rubber --unsafe` in the Makefile.

## Running the examples

The examples run with MATLAB (R2025b) and a worktree of the `accuracy` branch,
so that the main checkout of Dynare is left alone:

```sh
make dynare     # git worktree add ~/works/Dynare/accuracy accuracy, submodules, meson build
```

Variables (override on the command line): `DYNARE_REPO`
(`~/works/Dynare/dynare`), `DYNARE_ROOT` (`~/works/Dynare/accuracy`),
`DYNARE_BUILD_DIR` (`$(DYNARE_ROOT)/build-matlab`), `MATLAB_ROOT`
(`~/matlab/R2025b`). The branch changes neither the preprocessor nor the mex
files, so `DYNARE_BUILD_DIR` can point to the build of another checkout.

```sh
make results    # every example, about forty minutes
make growth     # one example: growth linear burnside dhm ep_burnside irr bs
make clean      # LaTeX and Dynare temporary files
make distclean  # also the PDF, the stamps and the .mat files of the results
```

A stamp in `mod/results/` records that an example ran with its current sources;
the stamps are not versioned, so after a checkout each example runs once. After
a run, commit the updated `mod/results/` (`.tex`, `.csv`, `dynare.hash`).

Indicative timings: a few minutes for `growth`, `linear`, `burnside`, `dhm` and
`ep_burnside`; about eight minutes for `irr`; about twenty-five minutes for `bs`,
most of it in the stochastic extended path.

- **Irreversible investment** (`irr`) chains its runs in one MATLAB session:
  `occbin_irr`, `ep_irr`, `sep_irr`, `sep_irr_h2` and `sep_irr_h4` (stochastic
  extended path of order 1, without and with the hybrid correction of order 2
  and 4) save snapshots of their solutions (`accuracy.context`), which
  `compare_irr` scores against the equations of the `min` formulation
  (`accuracy.compare`). The same comparison is repeated for a range of
  volatilities and two degrees of risk aversion (`sweep_irr`), the saved
  contexts being recalibrated in place (`recalibrate.m`).
- **Borrowing constraint** (`bs`) follows the same pattern (`occbin_bs`,
  `ep_bs`, `sep_bs`, `sep_bs_h2`, `sep_bs_h4`, `sep_bs_o2`, `sep_bs_o3`,
  `compare_bs`), with a global solution computed by the endogenous grid method
  (`buffer_egm.m`) as the reference. The rules are scored at states drawn from
  the ergodic distribution of that solution, with `score_rules.m`.

`compare_irr` uses the `log_variables` option of the Newton normalization, so it
needs a version of the branch that has it.

## Exploratory: New Keynesian model with a zero lower bound

`mod/zlb_*.inc` and `mod/occbin_zlb.mod` hold a New Keynesian model with
Epstein–Zin preferences and a zero lower bound on the nominal rate, streamlined
from `tests/ep/hybrid/rs.mod` of Dynare. It is not used by the slides and has no
Makefile target. The same economy is written three ways: with a `max`
(`zlb_max_model.inc`, the reference equations), as a complementarity condition
(`zlb_mcp_model.inc`, for the extended path solvers `lmmcp`,
`options_.ep.solve_algo = 10`, or PATH, `11`), and with two OccBin regimes
(`occbin_zlb.mod`). The calibration is set by macro variables in
`zlb_calibration.inc`.

PATH is not distributed with Dynare: the MATLAB interface is available at
<https://pages.cs.wisc.edu/~ferris/path.html>, its folder must be named
`mcppath` and be on the MATLAB path, and problems beyond 300 variables need the
licence string of <http://www.cs.wisc.edu/~ferris/path/LICENSE> in the
environment variable `PATH_LICENSE_STRING`. PATH reads its options from
`path.opt` in the current folder; `crash_method none` with
`nms_searchtype line` worked on this model.

## Licence

The code (`mod/`, `Makefile`, `.github/`) is under the GNU General Public License,
version 3 or later (`COPYING`); the slides (`tex/`) and the results they read are
dedicated to the public domain under CC0 1.0 (`tex/LICENSE`). See `LICENSE`.
