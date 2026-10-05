// The New Keynesian model with the zero lower bound written as a max: the
// yardstick of the comparison, whose residuals vanish at the bound. Newton fails
// at the kink of the max, so the extended paths solve zlb_mcp_model.inc instead.

@#include "zlb_declarations.inc"

model;
@#include "zlb_core.inc"
  [name='zlb']
  Int = max(Ishadow, 0);
end;

@#include "zlb_calibration.inc"
