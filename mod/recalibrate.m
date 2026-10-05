function ctx = recalibrate(ctx, names, values, Sigma_e)

% Changes parameters and the covariance of the innovations in a snapshot of a run.
%
% INPUTS
% - ctx         [struct]    context, see accuracy.context
% - names       [cell]      parameter names
% - values      [double]    their new values
% - Sigma_e     [double]    new covariance matrix of the innovations, [] to keep it
%
% OUTPUTS
% - ctx         [struct]    updated context
%
% REMARKS
% The parameters changed must leave the deterministic steady state unchanged.
% The first order decision rule stored in the context, which the OccBin adapter
% uses and which does not depend on Sigma_e, is recomputed; the extended path
% adapters rebuild everything else from M_ when they are created.

for i=1:length(names)
    ctx.M_.params(strcmp(ctx.M_.param_names, names{i})) = values(i);
end
if ~isempty(Sigma_e)
    ctx.M_.Sigma_e = Sigma_e;
end

opts = ctx.options_;
opts.order = 1;
opts.k_order_solver = false;
[dr, info] = resol(0, ctx.M_, opts, ctx.oo_.dr, ctx.oo_.steady_state, ...
                   ctx.oo_.exo_steady_state, ctx.oo_.exo_det_steady_state);
if info(1)
    error('recalibrate: no first order solution (info = %u).', info(1))
end
ctx.oo_.dr = dr;
