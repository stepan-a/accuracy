/* Irreversible investment, stochastic extended path of order 1, with the hybrid
   correction of order 4. */

@#define hybrid = 4
@#include "sep_irr.inc"

ctx_sep_h4 = accuracy.context(M_, options_, oo_, 'extended_path');
save('ctx_sep_h4.mat', 'ctx_sep_h4');
