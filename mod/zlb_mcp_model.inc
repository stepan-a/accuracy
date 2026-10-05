// The New Keynesian model with the zero lower bound written as a complementarity
// condition, for the mixed complementarity solvers of the extended path
// (options_.ep.solve_algo = 10, lmmcp, or 11, PATH). The residual of the zlb
// equation does not vanish at the bound, so the comparisons use the equations of
// zlb_max_model.inc.

@#include "zlb_declarations.inc"

model;
@#include "zlb_core.inc"
  [name='zlb']
  Int = Ishadow ⟂ Int > 0;
end;

@#include "zlb_calibration.inc"
