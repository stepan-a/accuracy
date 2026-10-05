/* Stochastic growth model: the residuals by quadrature, the evaluation points,
   and the normalizations. Adapted from tests/accuracy/steadystate.mod of the
   accuracy branch of Dynare. */

var c k a;
varexo e;

parameters alpha beta delta rho sigmae;

alpha  = 0.36;
beta   = 0.99;
delta  = 0.025;
rho    = 0.95;
sigmae = 0.02;

model;
  [name='euler', accuracy_units='c']
  1/c = beta*(1/c(+1))*(alpha*exp(a(+1))*k^(alpha-1) + 1 - delta);
  [name='budget', accuracy_units='k']
  c + k = exp(a)*k(-1)^alpha + (1-delta)*k(-1);
  [name='technology']
  a = rho*a(-1) + e;
end;

steady_state_model;
  a = 0;
  k = ((1/beta-1+delta)/alpha)^(1/(alpha-1));
  c = k^alpha - delta*k;
end;

shocks;
  var e = sigmae^2;
end;

steady;

stoch_simul(order=1, periods=0, irf=0, nomoments, noprint);

options_.nograph = true;
ic = find(strcmp(M_.endo_names, 'c'));
ieuler = get_equation_number_by_tag('euler', M_);
mac = struct();

/* 1. At the steady state: zero under certainty equivalence, a risk correction of
      order sigma^2 under quadrature. */

T = 20;
noshocks = zeros(T, M_.exo_nbr);
common = {'dhm', false, 'points', {'simulation'}, 'periods', T, 'drop', 0, 'shocks', noshocks};

oo_ = accuracy.check(M_, options_, oo_, common{:}, 'integration', 'Deterministic');
mac.growthDetResid = sci(max(abs(oo_.accuracy.simulation.residuals(:))));

oo_ = accuracy.check(M_, options_, oo_, common{:}, 'quadrature_nodes', 9);
r1 = max(abs(oo_.accuracy.simulation.residuals(:)));
iu = oo_.accuracy.simulation.units.equations==ieuler;
u1 = max(abs(oo_.accuracy.simulation.units.values(iu,:)));
n1 = max(abs(oo_.accuracy.simulation.newton.values(ic,:)));

Sigma_e = M_.Sigma_e;
M_.Sigma_e = Sigma_e/100;
oo_ = accuracy.check(M_, options_, oo_, common{:}, 'quadrature_nodes', 9);
r2 = max(abs(oo_.accuracy.simulation.residuals(:)));
M_.Sigma_e = Sigma_e;

mac.growthSSResid = sci(r1);

// The quadrature error: the residual at the steady state must not depend on the
// number of nodes once they are enough, the integrand being smooth.
qnames = {'Three', 'Five', 'Nine'};
qnodes = [3 5 9];
oo_ = accuracy.check(M_, options_, oo_, common{:}, 'quadrature_nodes', 21);
R21 = oo_.accuracy.simulation.residuals;
for i=1:length(qnodes)
    oo_ = accuracy.check(M_, options_, oo_, common{:}, 'quadrature_nodes', qnodes(i));
    mac.(['growthSSGap' qnames{i}]) = sci(max(abs(oo_.accuracy.simulation.residuals(:)-R21(:))));
end
mac.growthSSResidSmall = sci(r2);
mac.growthSSRatio = sprintf('%.1f', r1/r2);
mac.growthSSUnits = sci(100*u1, 3);
mac.growthSSNewton = sci(n1, 3);

// The three Jacobian modes of the Newton normalization.
jmodes = {'steadystate', 'mean', 'quadrature'};
nv = cell(3,1);
for i=1:3
    oo_ = accuracy.check(M_, options_, oo_, common{:}, 'quadrature_nodes', 9, ...
                         'normalization', {'raw', 'newton'}, 'newton_jacobian', jmodes{i});
    nv{i} = oo_.accuracy.simulation.newton.values;
end
mac.growthGapMean = sci(max(abs(nv{2}(:)-nv{1}(:)))/max(abs(nv{1}(:))));
mac.growthGapQuad = sci(max(abs(nv{3}(:)-nv{1}(:)))/max(abs(nv{1}(:))));

/* 2. The three kinds of evaluation points, with all the normalizations. */

oo_ = accuracy.check(M_, options_, oo_, 'dhm', false, 'points', {'simulation', 'grid', 'irf'}, ...
                     'periods', 1000, 'drop', 100, 'quadrature_nodes', 5, ...
                     'grid_points', 500, 'grid_pilot_periods', 5000, ...
                     'irf_periods', 40, 'irf_size', 3, 'seed', 1, 'newton_by_equation', true);
acc = oo_.accuracy;
save('results/growth_accuracy.mat', 'acc');

// Table of the residuals of the equations, simulation and grid points.
rows = {};
for j=1:M_.orig_endo_nbr
    rows{end+1} = {tex_name(acc.equations.names{j}), ...
                   lg(mean(abs(acc.simulation.residuals(j,:)))), lg(max(abs(acc.simulation.residuals(j,:)))), ...
                   lg(mean(abs(acc.grid.residuals(j,:)))), lg(max(abs(acc.grid.residuals(j,:))))};
end
write_tabular('results/growth_residuals.tex', 'lrrrr', ...
              {'Equation', 'mean (sim.)', 'max (sim.)', 'mean (grid)', 'max (grid)'}, rows);

// The same table after each normalization. Newton: one row per variable, in
// percent of the steady state (levels for a variable whose steady state is zero).
rows = {};
for j=1:M_.orig_endo_nbr
    ns = abs(acc.simulation.newton.values(j,:));
    ng = abs(acc.grid.newton.values(j,:));
    unit = '\%';
    if acc.simulation.newton.is_absolute(j)
        unit = 'level';
    end
    rows{end+1} = {sprintf('$%s$', M_.endo_names{j}), unit, lg(mean(ns)), lg(max(ns)), lg(mean(ng)), lg(max(ng))};
end
write_tabular('results/growth_newton_log.tex', 'llrrrr', ...
              {'Variable', 'unit', 'mean (sim.)', 'max (sim.)', 'mean (grid)', 'max (grid)'}, rows);

// Decomposition of the Newton step by equation, as returned by the Newton
// normalization with the newton_by_equation option: element (v,e,p) is the
// deviation of variable v the residual of equation e alone induces at point p,
// in percent of the steady state (levels for a). J0 is the Jacobian at the
// steady state, the default.
nw = acc.simulation.newton;
ys = oo_.steady_state;
neq = size(nw.by_equation, 2);
rows = {};
for e=1:neq
    dye = abs(squeeze(nw.by_equation(:,e,:)));
    if strcmp(acc.equations.names{e}, 'budget')
        mac.growthBudgetNewtonK = sprintf('%.2f', log10(mean(dye(strcmp(M_.endo_names, 'k'),:))));
    end
    row = {tex_name(acc.equations.names{e})};
    for v=1:M_.endo_nbr
        row = [row, {lg(mean(dye(v,:))), lg(max(dye(v,:)))}];
    end
    rows{end+1} = row;
end
dyt = abs(nw.values);
row = {'all (Newton step)'};
for v=1:M_.endo_nbr
    row = [row, {lg(mean(dyt(v,:))), lg(max(dyt(v,:)))}];
end
rows{end+1} = row;
gap = max(max(abs(squeeze(sum(nw.by_equation, 2))-nw.values)));
dprintf('growth: the decomposition by equation adds up to the Newton step up to %.2e.', gap);
header = {['& ' strjoin(cellfun(@(v) sprintf('\\multicolumn{2}{c}{$%s$}', v), M_.endo_names', ...
                                'UniformOutput', false), ' & ') ' \\ ' ...
           '\cmidrule(lr){2-3}\cmidrule(lr){4-5}\cmidrule(lr){6-7} Residual of']};
for v=1:M_.endo_nbr
    header = [header, {'mean', 'max'}];
end
write_tabular('results/growth_newton_by_equation.tex', 'lrrrrrr', header, rows);
mac.growthNewtonGap = sci(gap);
mac.growthSSc = sprintf('%.2f', ys(strcmp(M_.endo_names, 'c')));
mac.growthSSk = sprintf('%.0f', ys(strcmp(M_.endo_names, 'k')));

// Units of a variable: one row per tagged equation, in percent of the variable.
rows = {};
for j=1:length(acc.simulation.units.equations)
    us = 100*abs(acc.simulation.units.values(j,:));
    ug = 100*abs(acc.grid.units.values(j,:));
    rows{end+1} = {sprintf('%s [$%s$]', tex_name(acc.equations.names{acc.simulation.units.equations(j)}), ...
                           acc.simulation.units.names{j}), ...
                   lg(mean(us)), lg(max(us)), lg(mean(ug)), lg(max(ug))};
end
ib = acc.simulation.units.equations==get_equation_number_by_tag('budget', M_);
mac.growthBudgetUnitsK = sprintf('%.2f', log10(mean(100*abs(acc.simulation.units.values(ib,:)))));
write_tabular('results/growth_units_log.tex', 'lrrrr', ...
              {'Equation [variable]', 'mean (sim.)', 'max (sim.)', 'mean (grid)', 'max (grid)'}, rows);

// Table of the implied errors on the decision rule, Newton normalization, and
// the Euler equation error in units of consumption.
rows = {};
for j=1:M_.orig_endo_nbr
    nw = acc.simulation.newton;
    unit = '\%';
    if nw.is_absolute(j)
        unit = 'level';
    end
    rows{end+1} = {sprintf('$%s$', M_.endo_names{j}), unit, ...
                   sprintf('$%s$', sci(mean(abs(nw.values(j,:))))), sprintf('$%s$', sci(max(abs(nw.values(j,:)))))};
end
write_tabular('results/growth_newton.tex', 'llrr', {'Variable', 'unit', 'mean', 'max'}, rows);

un = acc.simulation.units;
mac.growthUnitsMean = sci(100*mean(abs(un.values(1,:))));
mac.growthUnitsMax = sci(100*max(abs(un.values(1,:))));
mac.growthNewtonCMean = sci(mean(abs(acc.simulation.newton.values(ic,:))));
mac.growthNewtonCMax = sci(max(abs(acc.simulation.newton.values(ic,:))));
mac.growthUnitsMeanLog = sprintf('%.2f', log10(100*mean(abs(un.values(1,:)))));
mac.growthUnitsMaxLog = sprintf('%.2f', log10(100*max(abs(un.values(1,:)))));
mac.growthNewtonCMeanLog = sprintf('%.2f', log10(mean(abs(acc.simulation.newton.values(ic,:)))));
mac.growthNewtonCMaxLog = sprintf('%.2f', log10(max(abs(acc.simulation.newton.values(ic,:)))));
mac.growthSimPoints = sprintf('%u', acc.simulation.npoints);
mac.growthGridPoints = sprintf('%u', acc.grid.npoints);
mac.growthIrfPoints = sprintf('%u', acc.irf.npoints);

// Grid map: lagged states and the Euler equation error in units of consumption.
sidx = oo_.dr.order_var(M_.nstatic+(1:M_.nspred));
S = acc.grid.points.state_lag;
ik = find(strcmp(M_.endo_names(sidx), 'k'));
ia = find(strcmp(M_.endo_names(sidx), 'a'));
write_csv('results/growth_grid.csv', {'k', 'a', 'lunits', 'lresid'}, ...
          [S(ik,:)' S(ia,:)' log10(abs(acc.grid.units.values(1,:))'+realmin) ...
           log10(abs(acc.grid.residuals(ieuler,:))'+realmin)]);
mac.growthGridKlo = sprintf('%.2f', acc.grid.points.bounds(ik,1));
mac.growthGridKhi = sprintf('%.2f', acc.grid.points.bounds(ik,2));
mac.growthGridAlo = sprintf('%.3f', acc.grid.points.bounds(ia,1));
mac.growthGridAhi = sprintf('%.3f', acc.grid.points.bounds(ia,2));

// Profile of the Euler equation error along the impulse response.
write_csv('results/growth_irf.csv', {'t', 'lunits', 'lresid'}, ...
          [acc.irf.points.irf_period(:) log10(abs(acc.irf.units.values(1,:))'+realmin) ...
           log10(abs(acc.irf.residuals(ieuler,:))'+realmin)]);

/* 3. The Euler equation error, in units of consumption, against the order of
      approximation, at the same simulated points. */

rows = {};
@#for order in 1:3
stoch_simul(order=@{order}, periods=0, irf=0, nomoments, noprint);
oo_ = accuracy.check(M_, options_, oo_, 'dhm', false, 'points', {'simulation'}, ...
                     'periods', 1000, 'drop', 100, 'quadrature_nodes', 5, 'seed', 1, 'pruning', false);
iu = oo_.accuracy.simulation.units.equations==ieuler;
u = 100*abs(oo_.accuracy.simulation.units.values(iu,:));
rows{end+1} = {'@{order}', lg(mean(u)), lg(max(u))};
@#endfor
write_tabular('results/growth_orders.tex', 'crr', ...
              {'Order', 'mean', 'max'}, rows);

write_macros('results/growth_macros.tex', mac);
