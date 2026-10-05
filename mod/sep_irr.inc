// Irreversible investment, stochastic extended path of order 1: the tree of
// future histories is built on a three node Gauss-Hermite rule. With a hybrid
// order h>1, the paths at the leaves of the tree are corrected by the constant of
// the perturbation approximation of order h. Included by sep_irr.mod,
// sep_irr_h2.mod and sep_irr_h4.mod, which set the macro variable hybrid.

@#include "irr_min_model.inc"

@#if hybrid > 2
// The order 4 correction is the constant g_0 of the k-order solver; the
// derivatives of the model are only computed up to the largest order requested.
stoch_simul(order=@{hybrid}, k_order_solver, periods=0, irf=0, nomoments, noprint);
stoch_simul(order=1, periods=0, irf=0, nomoments, noprint);
options_.k_order_solver = true;
@#endif
@#if hybrid == 2
stoch_simul(order=2, periods=0, irf=0, nomoments, noprint);
stoch_simul(order=1, periods=0, irf=0, nomoments, noprint);
@#endif

options_.ep.periods = 60;
options_.ep.stochastic.order = 1;
options_.ep.stochastic.IntegrationAlgorithm = 'Tensor-Gaussian-Quadrature';
options_.ep.stochastic.quadrature.nodes = 3;
options_.ep.stochastic.hybrid_order = @{hybrid};
