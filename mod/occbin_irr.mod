/* Irreversible investment, solved with OccBin. Derived from
   tests/occbin/model_irrcap_twoconstraints/dynrbc.mod, keeping only the
   irreversibility constraint, which is a genuine complementarity. The second
   constraint of that model is an asymmetric adjustment cost, a regime indicator
   rather than a complementarity, and has no counterpart in a min formulation.

   Copied from tests/accuracy/occbin_irr.mod. The same economy is written with a
   min in ep_irr.mod and sep_irr.mod. */

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

  [name='investment', bind='IRR']
  i - log(PHI*ziss) = 0;
  [name='investment', relax='IRR']
  lambdak = 0;

  [name='technology']
  a = RHOA*a(-1) + erra;
end;

occbin_constraints;
  name 'IRR'; bind i-steady_state(i)<log(PHI); relax lambdak<0;
end;

@#include "irr_calibration.inc"

options_.nograph = true;

stoch_simul(order=1, periods=0, irf=0, nomoments, noprint);

ctx_occbin = accuracy.context(M_, options_, oo_, 'occbin');
save('ctx_occbin.mat', 'ctx_occbin');
