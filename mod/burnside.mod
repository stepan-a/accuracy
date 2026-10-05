/* The model of Burnside (1998), whose exact solution is known in closed form.
   Adapted from tests/accuracy/burnside.mod. */

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

T = 200;
set_dynare_seed(1);
ex_ = 0.0348*randn(T, 1);

ieuler = get_equation_number_by_tag('euler', M_);
exact = @(S, U) burnside_exact_policy(S, U, M_, oo_);
common = {'dhm', false, 'points', {'simulation'}, 'periods', T, 'drop', 0, 'shocks', ex_};

// The exact solution: its residual is pure quadrature error.
nnodes = [3 5 7 9 11 13 15 21];
qerr = NaN(length(nnodes), 1);
for i=1:length(nnodes)
    oo_ = accuracy.check(M_, options_, oo_, common{:}, 'policy', exact, 'quadrature_nodes', nnodes(i));
    qerr(i) = max(abs(oo_.accuracy.simulation.residuals(ieuler,:)));
end

// The perturbation solutions: approximation error, measured with 21 nodes.
perr = NaN(4,1);
perr21 = NaN(4,1);
@#for order in 1:4
stoch_simul(order=@{order}, periods=0, irf=0, nomoments, noprint, k_order_solver);
oo_ = accuracy.check(M_, options_, oo_, common{:}, 'quadrature_nodes', 21);
perr(@{order}) = max(abs(oo_.accuracy.simulation.residuals(ieuler,:)));
nerr(@{order}) = max(abs(oo_.accuracy.simulation.newton.values(strcmp(M_.endo_names, 'y'),:)));
@#endfor

write_csv('results/burnside_nodes.csv', {'nodes', 'lerr'}, [nnodes(:) log10(qerr)]);
write_csv('results/burnside_orders.csv', {'order', 'lerr'}, [(1:4)' log10(perr)]);

mac = struct();
mac.burnsideExactFive = sci(qerr(nnodes==5));
mac.burnsideExactEleven = sci(qerr(nnodes==11));
mac.burnsideExactTwentyOne = sci(qerr(nnodes==21));
mac.burnsideOrderOne = sci(perr(1));
mac.burnsideOrderTwo = sci(perr(2));
mac.burnsideOrderThree = sci(perr(3));
mac.burnsideOrderFour = sci(perr(4));
mac.burnsideNewtonOne = sprintf('%.1f', nerr(1));
mac.burnsideNewtonFour = sprintf('%.1f', nerr(4));
mac.burnsideDecades = sprintf('%.1f', log10(perr(4)/qerr(end)));
mac.burnsideSteadyY = sprintf('%.2f', oo_.steady_state(strcmp(M_.endo_names, 'y')));
write_macros('results/burnside_macros.tex', mac);
