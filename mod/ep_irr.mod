/* Irreversible investment, deterministic extended path (stochastic order 0). */

@#include "irr_min_model.inc"

options_.ep.periods = 60;
options_.ep.stochastic.order = 0;

ctx_ep = accuracy.context(M_, options_, oo_, 'extended_path');
save('ctx_ep.mat', 'ctx_ep');
