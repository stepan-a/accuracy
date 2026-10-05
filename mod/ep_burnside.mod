/* The extended path, deterministic and stochastic, on the model of Burnside
   (1998). Adapted from tests/accuracy/ep_burnside.mod.

   Two quadratures are at play, and they are independent: the one the extended
   path uses to build its tree of future histories (options_.ep.stochastic), and
   the one the accuracy check uses to take the expectation of the residuals (the
   quadrature_nodes option). */

var y x;
varexo e;

parameters beta theta rho xbar;

xbar  =  0.0179;
rho   = -0.139;
theta = -10;
beta  =  0.95;

model;
  [name='euler']
  y = beta*exp(theta*x(+1))*(1+y(+1));
  [name='technology']
  x = (1-rho)*xbar + rho*x(-1) + e;
end;

shocks;
  var e; stderr 0.0348;
end;

initval;
  x = xbar;
  y = beta*exp(theta*xbar)/(1-beta*exp(theta*xbar));
end;

steady;

stoch_simul(order=1, periods=0, irf=0, nomoments, noprint);
options_.nograph = true;

ieuler = get_equation_number_by_tag('euler', M_);

T = 20;
set_dynare_seed(1);
ex_ = 0.0348*randn(T, 1);
common = {'points', {'simulation'}, 'periods', T, 'drop', 0, 'shocks', ex_, 'quadrature_nodes', 11, 'dhm', false};

// label, horizon, stochastic order, hybrid order. The hybrid correction adds, at
// the leaves of the tree, the constant of a perturbation approximation of order
// 2 (0.5*ghs2) or 4 (g_0): it only acts when the stochastic order is positive.
runs = {'EP', 40, 0, 0; ...
        'EP', 200, 0, 0; ...
        'SEP', 40, 1, 0; ...
        'SEP', 40, 1, 2; ...
        'SEP', 40, 1, 4; ...
        'SEP', 40, 2, 0; ...
        'SEP', 40, 2, 2; ...
        'SEP', 40, 2, 4; ...
        'SEP', 40, 3, 0};

rows = {};
err = NaN(size(runs,1), 1);
failures = zeros(size(runs,1), 1);
options_.ep.stochastic.IntegrationAlgorithm = 'Tensor-Gaussian-Quadrature';
options_.ep.stochastic.quadrature.nodes = 3;
options_.k_order_solver = true; // the order 4 hybrid correction needs dr.g_0
for i=1:size(runs,1)
    options_.ep.periods = runs{i,2};
    options_.ep.stochastic.order = runs{i,3};
    options_.ep.stochastic.hybrid_order = runs{i,4};
    t0 = tic;
    oo_ = accuracy.check(M_, options_, oo_, 'solution', 'extended_path', common{:});
    elapsed = toc(t0);
    err(i) = max(abs(oo_.accuracy.simulation.residuals(ieuler,:)));
    hybrid = '--';
    if runs{i,4}
        hybrid = sprintf('%u', runs{i,4});
    end
    rows{end+1} = {runs{i,1}, sprintf('%u', runs{i,2}), sprintf('%u', runs{i,3}), hybrid, ...
                   sprintf('$%s$', sci(err(i))), sprintf('%.1f', elapsed)};
    failures(i) = oo_.accuracy.simulation.eval_failures;
end
options_.ep.stochastic.hybrid_order = 0;
options_.k_order_solver = false;

// The same points, scored against the perturbation solutions and the exact one.
exact = @(S, U) burnside_exact_policy(S, U, M_, oo_);
oo_ = accuracy.check(M_, options_, oo_, 'policy', exact, common{:});
rows{end+1} = {'Exact', '--', '--', '--', ...
               sprintf('$%s$', sci(max(abs(oo_.accuracy.simulation.residuals(ieuler,:))))), '--'};

prows = {};
pmax = NaN(4,1);
@#for order in [1, 2, 4]
stoch_simul(order=@{order}, periods=0, irf=0, nomoments, noprint, k_order_solver);
oo_ = accuracy.check(M_, options_, oo_, common{:});
prows{end+1} = {'Perturbation @{order}', '--', '--', '--', ...
                sprintf('$%s$', sci(max(abs(oo_.accuracy.simulation.residuals(ieuler,:))))), '--'};
pmax(@{order}) = max(abs(oo_.accuracy.simulation.residuals(ieuler,:)));
@#endfor
rows = [prows, rows];

write_tabular('results/ep_burnside.tex', 'lcccrr', ...
              {'Solution', 'horizon', 'order', 'hybrid', 'max $|R|$', 'time (s)'}, rows);

mac = struct();
mac.epBurnsidePoints = sprintf('%u', T);
mac.epbFailures = sprintf('%u', sum(failures));
mac.epbDet = sci(err(1));
mac.epbSEPone = sci(err(3));
mac.epbSEPoneHtwo = sci(err(4));
mac.epbSEPoneHfour = sci(err(5));
mac.epbSEPtwo = sci(err(6));
mac.epbSEPtwoHfour = sci(err(8));
mac.epbSEPthree = sci(err(9));
mac.epbPertFour = sci(pmax(4));
write_macros('results/ep_burnside_macros.tex', mac);
