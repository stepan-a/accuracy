function Y = buffer_exact_policy(S, U, M_, oo_, sol)

% The global solution of the borrowing constraint model, as a decision rule for
% the accuracy checks (policy option).
%
% INPUTS
% - S           [double]    ns×N lagged states (b, y), in deviation from the
%                           steady state, in DR order
% - U           [double]    1×N innovations
% - M_          [struct]    describing the model
% - oo_         [struct]    results
% - sol         [struct]    solution returned by buffer_egm
%
% OUTPUTS
% - Y           [double]    n×N current endogenous variables, in declaration order
%
% REMARKS
% Besides debt and consumption, the rule must return the expected consumption ec
% and the multiplier lb, both functions of the expectation of next period's
% consumption, which is computed with the Gauss-Hermite rule of the solution.
% The multiplier is zero wherever the limit does not bind.

par = sol.par;
name = @(v) strcmp(M_.endo_names, v);
sidx = oo_.dr.order_var(M_.nstatic+(1:M_.nspred));
lag = repmat(oo_.dr.ys(sidx), 1, size(S,2))+S;
bprev = lag(strcmp(M_.endo_names(sidx), 'b'),:);
ylag = lag(strcmp(M_.endo_names(sidx), 'y'),:);

logy = par.RHO*log(ylag)+U;
[b, c] = buffer_policy(sol, bprev, logy);

N = length(b);
nq = length(sol.e);
ly1 = repmat(par.RHO*logy(:), 1, nq)+par.SIGMA*repmat(transpose(sol.e), N, 1);
[~, c1] = buffer_policy(sol, repmat(b, nq, 1), ly1(:));
c1 = reshape(c1, N, nq);

y = exp(logy(:));
lb = c.^(-par.GAMMAC)-par.BETA*par.R*(c1.^(-par.GAMMAC))*sol.w;
lb(b<par.M*y) = 0;

Y = zeros(M_.endo_nbr, N);
Y(name('b'),:) = b;
Y(name('c'),:) = c;
Y(name('ec'),:) = c1*sol.w;
Y(name('lb'),:) = lb;
Y(name('y'),:) = y;
