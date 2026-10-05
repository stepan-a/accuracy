/* New Keynesian model with Epstein-Zin preferences and a zero lower bound,
   solved with OccBin. The same economy is written with a max in
   zlb_max_model.inc and with a complementarity condition in zlb_mcp_model.inc. */

@#include "zlb_declarations.inc"

model;
@#include "zlb_core.inc"
  [name='zlb', relax='ZLB']
  Int = Ishadow;
  [name='zlb', bind='ZLB']
  Int = 0;
end;

occbin_constraints;
  name 'ZLB'; bind Ishadow<0; relax Ishadow>0;
end;

@#include "zlb_calibration.inc"

options_.nograph = true;

stoch_simul(order=1, periods=0, irf=0, nomoments, noprint);

ctx = accuracy.context(M_, options_, oo_, 'occbin');
save('ctx_occbin_zlb.mat', 'ctx');
