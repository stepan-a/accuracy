/* Irreversible investment, stochastic extended path of order 1. */

@#define hybrid = 0
@#include "sep_irr.inc"

ctx_sep = accuracy.context(M_, options_, oo_, 'extended_path');
save('ctx_sep.mat', 'ctx_sep');
