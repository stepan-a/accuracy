/* A linear model: its first order solution is exact, so every residual must
   vanish, at every kind of evaluation point. Adapted from tests/accuracy/linear.mod. */

var y pie r;
varexo e_y e_pie;

parameters delta sigma alpha kappa gamma1 gamma2;

delta  =  0.44;
kappa  =  0.18;
alpha  =  0.48;
sigma  = -0.06;
gamma1 =  1.5;
gamma2 =  0.5;

model;
  [name='IS curve']
  y = delta*y(-1) + (1-delta)*y(+1) + sigma*(r - pie(+1)) + e_y;
  [name='Phillips curve']
  pie = alpha*pie(-1) + (1-alpha)*pie(+1) + kappa*y + e_pie;
  [name='Taylor rule']
  r = gamma1*pie + gamma2*y;
end;

shocks;
  var e_y; stderr 0.63;
  var e_pie; stderr 0.4;
end;

steady;

stoch_simul(order=1, periods=0, irf=0, nomoments, noprint);
options_.nograph = true;

oo_ = accuracy.check(M_, options_, oo_, 'dhm', false, 'points', {'simulation', 'grid', 'irf'}, ...
                     'periods', 200, 'drop', 0, 'quadrature_nodes', 5, ...
                     'grid_points', 200, 'grid_pilot_periods', 2000, 'irf_periods', 20, 'seed', 1);

mac = struct();
names = {'Sim', 'Grid', 'Irf'};
kinds = {'simulation', 'grid', 'irf'};
for i=1:3
    mac.(['linearMax' names{i}]) = sci(max(abs(oo_.accuracy.(kinds{i}).residuals(:))));
    mac.(['linearN' names{i}]) = sprintf('%u', oo_.accuracy.(kinds{i}).npoints);
end
write_macros('results/linear_macros.tex', mac);
