/* Borrowing constraint model: deterministic and stochastic extended paths, and
   OccBin, scored against the min formulation at states drawn from the ergodic
   distribution of a global solution (endogenous grid method), which also
   provides the reference decision rule.

   The model below is a placeholder, Dynare needing something to process. */

@#include "bs_declarations.inc"

model;
  [name='expected consumption'] ec = c(1);
  [name='budget constraint'] c = y + b - R*b(-1);
  [name='borrowing'] lb = 0;
  [name='euler'] lb = 1/c^GAMMAC - BETA*R/c(+1)^GAMMAC;
  [name='income'] log(y) = RHO*log(y(-1)) + u;
end;

R = 1.05; BETA = 0.945; RHO = 0.9; SIGMA = 0.05; M = 1; GAMMAC = 1;

names = {'ep', 'sep', 'sep_h2', 'sep_h4', 'sep_o2', 'sep_o3', 'occbin'};
labels = {'EP', 'SEP', 'SEP, hybrid 2', 'SEP, hybrid 4', 'SEP, order 2', 'SEP, order 3', 'OccBin', 'Global (EGM)'};
tags = {'EP', 'SEP', 'SEPhTwo', 'SEPhFour', 'SEPoTwo', 'SEPoThree', 'OccBin', 'EGM'};
ctxs = cell(1, length(names));
for k=1:length(names)
    tmp = load(sprintf('ctx_%s_bs.mat', names{k}));
    ctxs{k} = tmp.ctx;
end
ref = ctxs{1};
Mr = ref.M_;
par = struct();
for pname = {'R', 'BETA', 'RHO', 'SIGMA', 'M', 'GAMMAC'}
    par.(pname{1}) = Mr.params(strcmp(Mr.param_names, pname{1}));
end

// Global solution.
t0 = tic;
// Built with many more Gauss-Hermite nodes than the checks use: a check with the
// solver's own rule would reproduce its expectation, and see only the
// interpolation error.
egm_nodes = 81;
sol = buffer_egm(par, struct('nq', egm_nodes));
elapsed_egm = toc(t0);
exact = @(S, U) buffer_exact_policy(S, U, ref.M_, ref.oo_, sol);

// States drawn from the ergodic distribution of the global solution.
P = 200;
burn = 1000;
set_dynare_seed(1);
innov = par.SIGMA*randn(burn+P, 1);
bp = par.M;
ly = 0;
bpath = zeros(burn+P, 1);
lypath = zeros(burn+P, 1);
for t=1:burn+P
    ly = par.RHO*ly+innov(t);
    bp = buffer_policy(sol, bp, ly);
    bpath(t) = bp;
    lypath(t) = ly;
end
sidx = ref.oo_.dr.order_var(Mr.nstatic+(1:Mr.nspred));
snames = Mr.endo_names(sidx);
S = zeros(length(sidx), P);
S(strcmp(snames, 'b'),:) = bpath(burn:burn+P-1);
S(strcmp(snames, 'y'),:) = exp(lypath(burn:burn+P-1));
x = transpose(innov(burn+1:burn+P));

t0 = tic;
scoring = {'ep_tolerance', 1e-12, 'occbin_periods', 100, 'occbin_maxit', 50};
out = score_rules(ref, [ctxs, {exact}], S, x, 'quadrature_nodes', 21, scoring{:});
elapsed = toc(t0);
save('results/compare_bs.mat', 'out', 'S', 'x');

K = length(labels);
ieq = @(n) find(strcmp(out.equations.names, n));
ieuler = ieq('euler');
iborrow = ieq('borrowing');
ibudget = ieq('budget constraint');
ib = strcmp(Mr.endo_names, 'b');
ic = strcmp(Mr.endo_names, 'c');
iy = strcmp(Mr.endo_names, 'y');
ybind = @(k) out.y{k}(ib,:)>=par.M*out.y{k}(iy,:)-1e-10;

fmean = @(v) mean(v(isfinite(v)));
fmax = @(v) max(v(isfinite(v)));
rows = {};
for k=1:K
    Rk = abs(out.residuals{k});
    db = abs(out.y{k}(ib,:)-out.y{K}(ib,:));
    dc = abs(out.y{k}(ic,:)-out.y{K}(ic,:));
    rows{end+1} = {labels{k}, lg(fmax(Rk(ieuler,:))), lg(fmean(Rk(ieuler,:))), lg(fmax(Rk(iborrow,:))), ...
                   sprintf('%.3f', fmean(db)), sprintf('%.3f', fmax(db)), sprintf('%.3f', fmean(dc)), ...
                   sprintf('%u', out.eval_failures(k))};
end
write_tabular('results/compare_bs.tex', 'lrrrrrrr', ...
              {['& \multicolumn{2}{c}{Euler, $\log_{10}|R|$} & & \multicolumn{2}{c}{$|b-b^\star|$} & & \\ ' ...
                '\cmidrule(lr){2-3}\cmidrule(lr){5-6} Solution'], 'max', 'mean', 'compl.', ...
               'mean', 'max', 'mean $|c-c^\star|$', 'failures'}, rows);

// The quadrature error of the check, on the global solution. With the solver's
// own rule (marked), the check reproduces the solution's expectation: zero where
// the limit binds, the multiplier being computed with that rule, and the
// interpolation error elsewhere. With more nodes, it converges to the true error
// of the solution; with fewer, it adds its own quadrature error.
qrows = {};
qnodes = [9 21 41 egm_nodes 161];
qmax = zeros(size(qnodes));
qmean = zeros(size(qnodes));
for q = 1:length(qnodes)
    o = score_rules(ref, {exact}, S, x, 'quadrature_nodes', qnodes(q));
    r = abs(o.residuals{1}(ieuler,:));
    qmax(q) = max(r);
    qmean(q) = mean(r);
    label = sprintf('%u', qnodes(q));
    if qnodes(q)==egm_nodes
        label = [label '$^\star$'];
    end
    qrows{end+1} = {label, lg(qmax(q)), lg(qmean(q))};
end
write_tabular('results/compare_bs_quadrature.tex', 'crr', {'nodes', 'max', 'mean'}, qrows);

// Ergodic behaviour: each rule simulated with the same innovations, from the
// steady state.
Tsim = 600;
dropsim = 100;
set_dynare_seed(2);
innov_sim = par.SIGMA*randn(Tsim, 1);
slack = NaN(Tsim, K);
t0 = tic;
for k=1:K
    if k<K
        Mk = ctxs{k}.M_;
    else
        Mk = Mr;
    end
    ysim = accuracy.simulate(out.rules{k}, Mk, Tsim, [], innov_sim);
    yref = NaN(Mr.endo_nbr, Tsim+1);
    yref(out.maps{k},:) = ysim;
    slack(:,k) = transpose(par.M*yref(iy,2:end)-yref(ib,2:end));
end
elapsed_sim = toc(t0);
srows = {};
for k=1:K
    v = slack(dropsim+1:end,k);
    srows{end+1} = {labels{k}, sprintf('%.3f', fmean(v)), sprintf('%.3f', fmax(v)), ...
                    sprintf('%.0f\\%%', 100*mean(v(isfinite(v))<1e-8)), sprintf('%u', sum(~isfinite(v)))};
end
write_tabular('results/compare_bs_ergodic.tex', 'lrrrr', ...
              {'Solution', 'mean buffer', 'max buffer', 'at the limit', 'failures'}, srows);
write_csv('results/compare_bs_sim.csv', [{'t'}, tags], [(1:Tsim)' slack]);

// Policy functions: debt chosen against debt inherited, at average income and
// no innovation.
bgrid = linspace(-1.5, 1.0, 101);
nb = length(bgrid);
L = zeros(length(sidx), nb);
L(strcmp(snames, 'b'),:) = bgrid;
L(strcmp(snames, 'y'),:) = 1;
pol = zeros(nb, K);
for k=1:K
    p = out.rules{k};
    map = out.maps{k};
    own = zeros(length(p.state_index), nb);
    for j=1:length(p.state_index)
        own(j,:) = L(strcmp(snames, Mr.endo_names{map(p.state_index(j))}),:);
    end
    Yk = p.eval(bsxfun(@minus, own, p.steady_state(p.state_index)), zeros(1, nb), []);
    yk = NaN(Mr.endo_nbr, nb);
    yk(map,:) = Yk;
    pol(:,k) = yk(ib,:)';
end
write_csv('results/compare_bs_policy.csv', [{'bprev'}, tags], [bgrid' pol]);

// Along the ergodic sample: distance to the limit under the global solution, and
// the Euler residual of each rule.
write_csv('results/compare_bs_points.csv', [{'t', 'slack'}, strcat({'r'}, tags)], ...
          [(1:P)' (par.M*out.y{K}(iy,:)'-out.y{K}(ib,:)') ...
           cell2mat(cellfun(@(r) log10(abs(r(ieuler,:))'+realmin), transpose(out.residuals), 'UniformOutput', false))]);

mac = struct();
mac.bsPoints = sprintf('%u', P);
mac.bsBurn = sprintf('%u', burn);
mac.bsEgmIterations = sprintf('%u', sol.iterations);
mac.bsEgmNodes = sprintf('%u', egm_nodes);
// Floor of the check with the 21 nodes of the comparison, and gap between the
// mean Euler residuals of the decision rules and of the global solution, in
// orders of magnitude.
mac.bsQuadFloorMax = sprintf('%.1f', log10(qmax(qnodes==21)));
mac.bsQuadFloorMean = sprintf('%.1f', log10(qmean(qnodes==21)));
gap = arrayfun(@(k) log10(fmean(abs(out.residuals{k}(ieuler,:)))/fmean(abs(out.residuals{K}(ieuler,:)))), 1:K-1);
mac.bsGapMin = sprintf('%.1f', min(gap));
mac.bsGapMax = sprintf('%.1f', max(gap));
mac.bsEgmTime = sprintf('%.0f', elapsed_egm);
mac.bsElapsed = sprintf('%.0f', elapsed);
mac.bsElapsedSim = sprintf('%.0f', elapsed_sim);
mac.bsTsim = sprintf('%u', Tsim-dropsim);
mac.bsFailures = sprintf('%u', sum(out.eval_failures));
for k=1:K
    mac.(['bsEuler' tags{k}]) = sci(fmax(abs(out.residuals{k}(ieuler,:))));
    mac.(['bsEulerMean' tags{k}]) = sci(fmean(abs(out.residuals{k}(ieuler,:))));
    mac.(['bsBuffer' tags{k}]) = sprintf('%.2f', fmean(slack(dropsim+1:end,k)));
    mac.(['bsAtLimit' tags{k}]) = sprintf('%.0f', 100*mean(slack(dropsim+1:end,k)<1e-8));
end
write_macros('results/compare_bs_macros.tex', mac);
