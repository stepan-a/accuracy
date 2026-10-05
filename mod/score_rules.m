function out = score_rules(ref, rules, state_lag, x, varargin)

% Scores several decision rules at given states, against the equations of one
% context.
%
% INPUTS
% - ref         [struct]    context providing the equations, see accuracy.context
% - rules       [cell]      each element is either a context, or a function
%                           handle Y = f(S,U) standing for a decision rule of the
%                           reference model (see the policy option of accuracy.check)
% - state_lag   [double]    ns×P lagged states, in levels, in the DR order of the
%                           states of the reference model
% - x           [double]    q×P innovations
% - varargin                accuracy options
%
% OUTPUTS
% - out         [struct]    residuals, y (reference ordering), regime,
%                           eval_failures, rules (the adapters) and maps (their
%                           variables in the reference ordering), one cell per
%                           rule, plus equations; newton, one cell per rule,
%                           when 'newton' is among the normalizations (see
%                           accuracy.newton_normalization)
%
% REMARKS
% This is what accuracy.compare does, except that the states are provided rather
% than drawn from a simulation of the reference rule. It is needed when the
% reference rule is not one of the rules compared, here a global solution
% computed outside of Dynare, and when the rules compared have ergodic sets far
% from each other.

opt = accuracy.default_options(ref.M_, ref.options_);
opt.normalization = {'raw'};
opt = accuracy.set_options(opt, varargin{:});

rm = accuracy.residual_model(ref.M_, ref.oo_);
out.equations = accuracy.equation_labels(ref.M_);
[nodes, weights] = accuracy.integration_nodes(ref.M_, opt);

kopt = opt;
kopt.solution = ref.solution;
kopt.policy = [];
pref = accuracy.policy_rule(ref.M_, ref.options_, ref.oo_, kopt);

K = length(rules);
out.residuals = cell(K, 1);
out.y = cell(K, 1);
out.regime = cell(K, 1);
out.eval_failures = zeros(K, 1);
out.rules = cell(K, 1);
out.maps = cell(K, 1);

for k=1:K
    kopt = opt;
    if isa(rules{k}, 'function_handle')
        kopt.policy = rules{k};
        p = accuracy.policy_rule(ref.M_, ref.options_, ref.oo_, kopt);
        map = transpose(1:ref.M_.endo_nbr);
    else
        kopt.solution = rules{k}.solution;
        kopt.policy = [];
        p = accuracy.policy_rule(rules{k}.M_, rules{k}.options_, rules{k}.oo_, kopt);
        map = accuracy.map_context(rules{k}, ref);
    end
    [ok, where] = ismember(map(p.state_index), pref.state_index);
    if ~all(ok)
        error('score_rules: the models do not have the same state variables.')
    end
    S = state_lag(where,:);
    [Y, ~, REG] = p.eval(bsxfun(@minus, S, p.steady_state(p.state_index)), x, []);
    pts = struct('description', 'given states', 'npoints', size(S, 2), ...
                 'state_lag', S, 'y', Y, 'x', x, 'latent', [], 'regime', REG, ...
                 'innovations', [], 'path', [], 'path_latent', [], 'bounds', [], ...
                 'irf_shock', [], 'irf_period', [], 'kind', 'shared', 'state_index', p.state_index);
    rmk = rm;
    rmk.map = map;
    [R, info, DY] = accuracy.residuals(p, rmk, opt, pts, nodes, weights);
    out.residuals{k} = R;
    if ismember('newton', opt.normalization)
        out.newton{k} = accuracy.newton_normalization(rmk, opt, R, DY, info.regimes, info.newton_by_equation);
    end
    out.y{k} = NaN(ref.M_.endo_nbr, size(S, 2));
    out.y{k}(map,:) = Y;
    out.regime{k} = REG;
    out.eval_failures(k) = info.eval_failures;
    out.rules{k} = p;
    out.maps{k} = map;
end
