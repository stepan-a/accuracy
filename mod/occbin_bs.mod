/* Borrowing constraint model, solved with OccBin. The model and the constraint
   are those of tests/occbin/model_borrcon. */

@#include "bs_declarations.inc"

model;
  [name='expected consumption']
  ec = c(1);
  [name='budget constraint']
  c = y + b - R*b(-1);
  [name='borrowing', bind='borrcon']
  lb = 0;
  [name='borrowing', relax='borrcon']
  b = M*y;
  [name='euler']
  lb = 1/c^GAMMAC - BETA*R/c(+1)^GAMMAC;
  [name='income']
  log(y) = RHO*log(y(-1)) + u;
end;

occbin_constraints;
  name 'borrcon'; bind lb<0; relax b>M*y;
end;

@#include "bs_calibration.inc"

options_.nograph = true;
stoch_simul(order=1, periods=0, irf=0, nomoments, noprint);

ctx = accuracy.context(M_, options_, oo_, 'occbin');
save('ctx_occbin_bs.mat', 'ctx');
