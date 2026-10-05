/* Borrowing constraint model, deterministic extended path. */

@#include "bs_min_model.inc"

stoch_simul(order=1, periods=0, irf=0, nomoments, noprint);
options_.ep.stochastic.order = 0;

ctx = accuracy.context(M_, options_, oo_, 'extended_path');
save('ctx_ep_bs.mat', 'ctx');
