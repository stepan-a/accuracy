/* Deterministic extended path, stochastic extended path and OccBin on the
   irreversible investment economy, all scored against the equations of the min
   formulation. Adapted from tests/accuracy/compare_irr.mod.

   The model below is a placeholder, Dynare needing something to process: the
   solutions compared are the snapshots saved by occbin_irr.mod, ep_irr.mod and
   sep_irr.mod. */

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

load('ctx_ep.mat');
load('ctx_sep.mat');
load('ctx_sep_h2.mat');
load('ctx_sep_h4.mat');
load('ctx_occbin.mat');

labels = {'EP', 'SEP', 'SEP, hybrid 2', 'SEP, hybrid 4', 'OccBin'};
tags = {'EP', 'SEP', 'SEPhTwo', 'SEPhFour', 'OccBin'};
P = 200;
set_dynare_seed(1);
innovations = 0.015*randn(P+50, 1);

t0 = tic;
out = accuracy.compare({ctx_ep, ctx_sep, ctx_sep_h2, ctx_sep_h4, ctx_occbin}, ...
                       'points', {'simulation'}, 'periods', P, 'drop', 50, 'shocks', innovations, ...
                       'max_points', P, 'quadrature_nodes', 5, 'occbin_periods', 60, ...
                       'ep_tolerance', 1e-12, 'dhm', false, 'noprint', false, ...
                       'normalization', {'raw', 'newton'}, 'newton_by_equation', true, ...
                       'log_variables', {'c', 'i', 'k'});
elapsed = toc(t0);
save('results/compare_irr.mat', 'out');

names = out.equations.names;
K = length(labels);
iocc = K;
ieuler = find(strcmp(names, 'euler'));
ib = find(strcmp(names, 'budget constraint'));
iid = find(strcmp(names, 'investment definition'));
iinv = find(strcmp(names, 'investment'));

// One row per method: base ten logarithm of the largest and mean absolute
// residual of the Euler equation, of the largest residual of the equations
// without expectation, and of the complementarity.
rows = {};
for k=1:K
    R = abs(out.residuals{k});
    rows{end+1} = {labels{k}, lg(max(R(ieuler,:))), lg(mean(R(ieuler,:))), ...
                   lg(max(max(R([ib iid],:)))), lg(max(R(iinv,:)))};
end
write_tabular('results/compare_irr_residuals.tex', 'lrrrr', ...
              {['& \multicolumn{2}{c}{Euler equation} & & \\ \cmidrule(lr){2-3} Solution'], ...
               'max', 'mean', 'budget, inv.~def.', 'complementarity'}, rows);

// Largest distance of each rule to the deterministic extended path, at the same
// states. The variables are in logs, so the distances are relative deviations.
endo = ctx_ep.M_.endo_names;
rows = {};
for k=2:K
    row = {labels{k}};
    for v = {'c', 'i', 'k', 'lambdak'}
        row{end+1} = sprintf('$%s$', sci(out.policy_distance(strcmp(endo, v{1}),k)));
    end
    rows{end+1} = row;
end
write_tabular('results/compare_irr_distance.tex', 'lrrrr', ...
              {'Distance to EP', '$c$', '$i$', '$k$', '$\lambda_k$'}, rows);

// Point by point: investment chosen by each rule, in percent of its steady
// state, the regime of OccBin, and the Euler equation residual.
ii = strcmp(endo, 'i');
iss = ctx_ep.oo_.steady_state(ii);
reg = out.regime{iocc};
if iscell(reg)
    reg = cell2mat(reg);
end
pdata = [(1:out.npoints)' 100*(out.y{1}(ii,:)'-iss) 100*(out.y{iocc}(ii,:)'-iss) reg(:)];
for k=1:K
    pdata = [pdata log10(abs(out.residuals{k}(ieuler,:))'+realmin)];
end
write_csv('results/compare_irr_points.csv', ...
          {'t', 'iep', 'iocc', 'regime', 'rep', 'rsep', 'rseph2', 'rseph4', 'rocc'}, pdata);

// Newton step on the current variables, the Jacobian being evaluated at each
// point: log10 of the implied error, in percent of the level (c, i and k are
// logarithms, see the log_variables option).
rows = {};
for k=1:K
    row = {labels{k}};
    for v = {'c', 'i', 'k'}
        dv = abs(out.newton{k}.values(strcmp(endo, v{1}),:));
        row = [row, {lg(mean(dv)), lg(max(dv))}];
    end
    rows{end+1} = row;
end
write_tabular('results/compare_irr_newton.tex', 'lrrrrrr', ...
              {['& \multicolumn{2}{c}{$c$} & \multicolumn{2}{c}{$i$} & \multicolumn{2}{c}{$k$} \\ ' ...
                '\cmidrule(lr){2-3}\cmidrule(lr){4-5}\cmidrule(lr){6-7} Solution'], ...
               'mean', 'max', 'mean', 'max', 'mean', 'max'}, rows);

// The same, equation by equation, for EP and OccBin: the error each equation
// alone induces on c, i and k, log10 of the mean in percent of the level.
rows = {};
for k=[1 iocc]
    for e=1:length(names)
        % The technology equation and the complementarity hold exactly for every
        % method: their rows would be empty.
        if out.equations.is_auxiliary(e) || ismember(names{e}, {'technology', 'investment'})
            continue
        end
        row = {labels{k}, tex_name(names{e})};
        for v = {'c', 'i', 'k'}
            d = abs(squeeze(out.newton{k}.by_equation(strcmp(endo, v{1}), e, :)));
            row = [row, {lg(mean(d)), lg(max(d))}];
        end
        rows{end+1} = row;
    end
end
write_tabular('results/compare_irr_newton_by_equation.tex', 'llrrrrrr', ...
              {['& & \multicolumn{2}{c}{$c$} & \multicolumn{2}{c}{$i$} & \multicolumn{2}{c}{$k$} \\ ' ...
                '\cmidrule(lr){3-4}\cmidrule(lr){5-6}\cmidrule(lr){7-8} Solution'], 'Residual of', ...
               'mean', 'max', 'mean', 'max', 'mean', 'max'}, rows);

mac = struct();
mac.irrNewtonOccI = sprintf('%.2f', max(abs(out.newton{iocc}.values(ii,:))));
mac.irrNewtonEPI = sprintf('%.2f', max(abs(out.newton{1}.values(ii,:))));
ic = strcmp(endo, 'c');
occ_by = @(eqname) max(abs(squeeze(out.newton{iocc}.by_equation(ic, strcmp(names, eqname), :))));
mac.irrByEqOccCEuler = sprintf('%.2f', occ_by('euler'));
mac.irrByEqOccCBudget = sprintf('%.2f', occ_by('budget constraint'));
mac.irrByEqOccCInvdef = sprintf('%.2f', occ_by('investment definition'));
mac.irrByEqOccCTotal = sprintf('%.2f', max(abs(out.newton{iocc}.values(ic,:))));
// For OccBin, the weight of the budget constraint and of the definition of
// investment relative to the Euler equation, on c, i and k: range of the ratios
// of the mean and of the max errors over the points.
ratios = struct('mean', [], 'max', []);
ie = strcmp(names, 'euler');
for v = {'c', 'i', 'k'}
    iv = strcmp(endo, v{1});
    de = abs(squeeze(out.newton{iocc}.by_equation(iv, ie, :)));
    for eqname = {'budget constraint', 'investment definition'}
        d = abs(squeeze(out.newton{iocc}.by_equation(iv, strcmp(names, eqname{1}), :)));
        ratios.mean(end+1) = mean(d)/mean(de);
        ratios.max(end+1) = max(d)/max(de);
    end
end
mac.irrByEqMeanRatioMin = sprintf('%.0f', min(ratios.mean));
mac.irrByEqMeanRatioMax = sprintf('%.0f', max(ratios.mean));
mac.irrByEqMaxRatioMin = sprintf('%.0f', min(ratios.max));
mac.irrByEqMaxRatioMax = sprintf('%.0f', max(ratios.max));
// Gain of the hybrid correction on the mean errors of EP, in orders of magnitude,
// averaged over c, i and k.
gain = [];
for v = {'c', 'i', 'k'}
    iv = strcmp(endo, v{1});
    gain(end+1) = log10(mean(abs(out.newton{1}.values(iv,:)))/mean(abs(out.newton{4}.values(iv,:))));
end
mac.irrHybridGain = sprintf('%.1f', mean(gain));
mac.irrPoints = sprintf('%u', out.npoints);
mac.irrNodes = sprintf('%u', out.nodes.nnodes);
mac.irrBinding = sprintf('%u', sum(reg(:)~=0));
mac.irrFloor = sprintf('%.2f', 100*log(PHI));
mac.irrCrossCheck = sci(out.cross_check(iocc));
mac.irrFailures = sprintf('%u', sum(out.eval_failures));
mac.irrElapsed = sprintf('%.0f', elapsed);
for k=1:K
    mac.(['irrEuler' tags{k}]) = sci(max(abs(out.residuals{k}(ieuler,:))));
    mac.(['irrEulerMean' tags{k}]) = sci(mean(abs(out.residuals{k}(ieuler,:))));
end
meanE = @(k) mean(abs(out.residuals{k}(ieuler,:)));
maxE = @(k) max(abs(out.residuals{k}(ieuler,:)));
mac.irrRatioOccEP = sprintf('%.1f', maxE(iocc)/maxE(1));
mac.irrRatioMeanOccEP = sprintf('%.1f', meanE(iocc)/meanE(1));
mac.irrRatioEPSEP = sprintf('%.1f', meanE(1)/meanE(2));
mac.irrRatioEPSEPhTwo = sprintf('%.1f', meanE(1)/meanE(3));
mac.irrRatioEPSEPhFour = sprintf('%.1f', meanE(1)/meanE(4));
mac.irrRatioMaxEPSEPhFour = sprintf('%.1f', maxE(1)/maxE(4));
mac.irrDistOccI = sprintf('%.2f', 100*out.policy_distance(ii,iocc));
mac.irrDistSepI = sprintf('%.2f', 100*out.policy_distance(ii,2));
mac.irrDistSepHFourI = sprintf('%.2f', 100*out.policy_distance(ii,4));
write_macros('results/compare_irr_macros.tex', mac);
