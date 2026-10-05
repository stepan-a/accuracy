/* Den Haan and Marcet test on the model of Burnside (1998): exact solution
   against a first order approximation. Adapted from tests/accuracy/dhm.mod. */

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

exact = @(S, U) burnside_exact_policy(S, U, M_, oo_);
common = {'points', {}, 'dhm_periods', 1000, 'dhm_replic', 200, 'dhm_alpha', 0.05, 'seed', 1};

oo_ = accuracy.check(M_, options_, oo_, 'policy', exact, common{:});
dexact = oo_.accuracy.dhm;

oo_ = accuracy.check(M_, options_, oo_, common{:});
dorder1 = oo_.accuracy.dhm;

stoch_simul(order=2, periods=0, irf=0, nomoments, noprint);
oo_ = accuracy.check(M_, options_, oo_, common{:});
dorder2 = oo_.accuracy.dhm;

write_csv('results/dhm_J.csv', {'exact', 'order1', 'order2'}, ...
          [dexact.J(1,:)' dorder1.J(1,:)' dorder2.J(1,:)']);

mac = struct();
mac.dhmDf = sprintf('%u', dexact.df(1));
mac.dhmInstruments = dexact.instruments_description;
labels = {'Exact', 'One', 'Two'};
tests = {dexact, dorder1, dorder2};
for i=1:3
    r = tests{i}.replications;
    mac.(['dhmMean' labels{i}]) = sprintf('%.2f', r.mean_J(1));
    mac.(['dhmLower' labels{i}]) = sprintf('%.3f', r.rejection_frequency_lower(1));
    mac.(['dhmUpper' labels{i}]) = sprintf('%.3f', r.rejection_frequency_upper(1));
end
// Instruments against sample length, order 2: rejection frequency in the upper
// tail, monomials of the normalized state up to degree d, T periods.
Ts = [1000 10000];
degrees = 1:3;
rows = {};
for d = degrees
    row = {sprintf('%u', d), sprintf('%u', d+1)};
    for T = Ts
        oo_ = accuracy.check(M_, options_, oo_, 'points', {}, 'dhm_periods', T, 'dhm_replic', 200, ...
                             'dhm_alpha', 0.05, 'seed', 1, 'dhm_instrument_degree', d);
        r = oo_.accuracy.dhm.replications;
        row = [row, {sprintf('%.2f', r.mean_J(1)), sprintf('%.0f\\%%', 100*r.rejection_frequency_upper(1))}];
    end
    rows{end+1} = row;
end
write_tabular('results/dhm_instruments.tex', 'cccccc', ...
              {['& & \multicolumn{2}{c}{$T=1000$} & \multicolumn{2}{c}{$T=10000$} \\ ' ...
                '\cmidrule(lr){3-4}\cmidrule(lr){5-6} degree $d$'], '$m$', 'mean $J$', 'rejected', ...
               'mean $J$', 'rejected'}, rows);
mac.dhmInstrumentsReplic = '200';

write_macros('results/dhm_macros.tex', mac);
