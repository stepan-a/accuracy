/* Irreversible investment, stochastic extended path of order 1, with the hybrid
   correction of order 2. */

@#define hybrid = 2
@#include "sep_irr.inc"

ctx_sep_h2 = accuracy.context(M_, options_, oo_, 'extended_path');
save('ctx_sep_h2.mat', 'ctx_sep_h2');
