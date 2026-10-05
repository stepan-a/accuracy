/* Irreversible investment: accuracy of EP, SEP with hybrid correction and OccBin as the volatility of
   technology rises, for two degrees of risk aversion. The snapshots saved by
   occbin_irr.mod, ep_irr.mod and sep_irr*.mod are recalibrated in place: GAMMAC
   and the variance of the innovations leave the steady state unchanged. The model
   below is a placeholder. */

@#include "irr_declarations.inc"

model;
  [name='euler'] exp(c) = exp(c(-1));
  [name='budget constraint'] exp(k) = exp(k(-1));
  [name='investment definition'] exp(i) = exp(i(-1));
  [name='investment'] lambdak = 0;
  [name='technology'] a = RHOA*a(-1) + erra;
end;

BETA = 0.96; ALPHA = 0.33; DELTAK = 0.10; GAMMAC = 2; RHOA = 0.9; PHI = 0.975;

shocks;
  var erra; stderr 0.015;
end;

// SEP without correction is close to EP, and hybrid 4 to hybrid 2, on this
// model (see compare_irr.mod): they are left out, the stochastic extended path
// being much slower when the innovations are large.
base = {load('ctx_ep.mat'), load('ctx_sep_h2.mat'), load('ctx_occbin.mat')};
base = {base{1}.ctx_ep, base{2}.ctx_sep_h2, base{3}.ctx_occbin};
tags = {'ep', 'seph2', 'occ'};
K = length(base);

gammas = [2 5];
sigmas = [0.005 0.01 0.015 0.02 0.03];
P = 60;
set_dynare_seed(1);
z = randn(P+50, 1);

results = [];
t0 = tic;
for g = gammas
    for s = sigmas
        ctxs = cellfun(@(c) recalibrate(c, {'GAMMAC'}, g, s^2), base, 'UniformOutput', false);
        out = accuracy.compare(ctxs, 'points', {'simulation'}, 'periods', P, 'drop', 50, 'shocks', s*z, ...
                               'max_points', P, 'quadrature_nodes', 5, 'occbin_periods', 60, ...
                               'ep_tolerance', 1e-12, 'dhm', false, 'noprint', true);
        ieuler = find(strcmp(out.equations.names, 'euler'));
        reg = out.regime{K};
        if iscell(reg)
            reg = cell2mat(reg);
        end
        row = [g s mean(reg(:)~=0)];
        for k=1:K
            r = abs(out.residuals{k}(ieuler,:));
            row = [row log10(mean(r(isfinite(r)))) log10(max(r(isfinite(r)))) out.eval_failures(k)];
        end
        results = [results; row];
        dprintf('sweep_irr: gamma %g, sigma %.3f, binding %.2f, mean log10 |R| %s, failures %s', g, s, row(3), ...
                mat2str(row(4:3:end), 3), mat2str(row(6:3:end)));
    end
end
elapsed = toc(t0);

header = {'gamma', 'sigma', 'binding'};
for k=1:K
    header = [header, strcat(tags{k}, {'mean', 'max', 'fail'})];
end
write_csv('results/sweep_irr.csv', header, results);

mac = struct();
mac.sweepPoints = sprintf('%u', P);
mac.sweepElapsed = sprintf('%.0f', elapsed);
mac.sweepFailures = sprintf('%u', sum(sum(results(:,6:3:end))));
mac.sweepSigmaMax = sprintf('%.2f', max(sigmas));
// Average gaps, in decades of the mean Euler residual, over the volatilities.
names = {'Two', 'Five'};
for i=1:2
    rg = results(:,1)==gammas(i);
    mac.(['sweepOccGap' names{i}]) = sprintf('%.1f', mean(results(rg,10)-results(rg,4)));
    mac.(['sweepHybridGain' names{i}]) = sprintf('%.1f', mean(results(rg,4)-results(rg,7)));
end
write_macros('results/sweep_irr_macros.tex', mac);
