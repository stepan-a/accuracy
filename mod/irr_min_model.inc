// The irreversible investment economy, with the complementarity written as a
// single nonlinear equation: the formulation the extended path can solve, and
// the one that provides the yardstick of the comparison. Shared by ep_irr.mod
// and sep_irr.mod, taken from tests/accuracy/ep_irr.mod.

@#include "irr_declarations.inc"

model;
# zkss = ((1/BETA-1+DELTAK)/ALPHA)^(1/(ALPHA-1));
# ziss = DELTAK*zkss;

  [name='euler', accuracy_units='c']
  -exp(c)^(-GAMMAC) + BETA*exp(c(1))^(-GAMMAC)*(1-DELTAK+ALPHA*exp(a(1))*exp(k)^(ALPHA-1))
    = -lambdak + BETA*(1-DELTAK)*lambdak(1);

  [name='budget constraint']
  exp(c) + exp(k) - (1-DELTAK)*exp(k(-1)) = exp(a)*exp(k(-1))^ALPHA;

  [name='investment definition']
  exp(i) = exp(k) - (1-DELTAK)*exp(k(-1));

  [name='investment']
  min(lambdak, i - log(PHI*ziss)) = 0;

  [name='technology']
  a = RHOA*a(-1) + erra;
end;

@#include "irr_calibration.inc"

options_.nograph = true;

stoch_simul(order=1, periods=0, irf=0, nomoments, noprint);
