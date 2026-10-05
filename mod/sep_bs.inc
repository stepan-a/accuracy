@#ifndef sep_order
  @#define sep_order = 1
@#endif

// Borrowing constraint model, stochastic extended path (order sep_order, 1 by
// default) on a three node Gauss-Hermite rule, with the hybrid correction of
// order hybrid when it is larger than one. Included by sep_bs.mod, sep_bs_h2.mod,
// sep_bs_h4.mod, sep_bs_o2.mod and sep_bs_o3.mod.

@#include "bs_min_model.inc"

@#if hybrid > 1
// The derivatives of the model are only computed up to the largest order requested.
stoch_simul(order=@{hybrid}, k_order_solver, periods=0, irf=0, nomoments, noprint);
@#endif
stoch_simul(order=1, periods=0, irf=0, nomoments, noprint);
options_.k_order_solver = @{hybrid} > 2;

options_.ep.stochastic.order = @{sep_order};
options_.ep.stochastic.IntegrationAlgorithm = 'Tensor-Gaussian-Quadrature';
options_.ep.stochastic.quadrature.nodes = 3;
options_.ep.stochastic.hybrid_order = @{hybrid};

ctx = accuracy.context(M_, options_, oo_, 'extended_path');
save('ctx_sep@{suffix}_bs.mat', 'ctx');
