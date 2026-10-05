function sol = buffer_egm(par, opt)

% Global solution of the borrowing constraint model by the endogenous grid method.
%
% INPUTS
% - par         [struct]    R, BETA, RHO, SIGMA, M, GAMMAC
% - opt         [struct]    optional: ny (nodes for log income), nk (nodes for
%                           the end of period debt), nq (Gauss-Hermite nodes),
%                           width (half width of the income grid, in unconditional
%                           standard deviations), blo (lowest debt), tol, maxit
%
% OUTPUTS
% - sol         [struct]    the solution, to be evaluated by buffer_policy
%
% REMARKS
% The model is
%     c(t) = y(t) + b(t) - R b(t-1)
%     c(t)^(-GAMMAC) = BETA R E_t[c(t+1)^(-GAMMAC)] + lb(t)
%     b(t) <= M y(t),  lb(t) >= 0,  lb(t) (M y(t) - b(t)) = 0
%     log y(t) = RHO log y(t-1) + u(t),  u ~ N(0, SIGMA^2)
% where b is debt. For each node of log income, the endogenous grid method maps a
% grid of end of period debt b(t), whose top point is the borrowing limit, onto
% the debt b(t-1) that makes the Euler equation hold with equality. Above the
% last of these points the constraint binds. The policy at an arbitrary income
% is obtained by linear interpolation in log income of the unconstrained policy,
% capped at the borrowing limit, so that the complementarity holds exactly.

if nargin<2
    opt = struct();
end
def = struct('ny', 101, 'nk', 600, 'nq', 21, 'width', 5, 'blo', -6, 'tol', 1e-13, 'maxit', 5000);
names = fieldnames(def);
for i=1:length(names)
    if ~isfield(opt, names{i})
        opt.(names{i}) = def.(names{i});
    end
end

sy = par.SIGMA/sqrt(1-par.RHO^2);
sol.par = par;
sol.logy = linspace(-opt.width*sy, opt.width*sy, opt.ny);
[e, w] = gauss_hermite(opt.nq);
sol.e = e;
sol.w = w;

% Grid of end of period debt, one per income node, denser near the limit.
s = linspace(0, 1, opt.nk);
bnext = zeros(opt.nk, opt.ny);
for i=1:opt.ny
    top = par.M*exp(sol.logy(i));
    bnext(:,i) = opt.blo+(top-opt.blo)*(1-(1-s').^2);
end
sol.bnext = bnext;

% Initial guess: debt at the limit, whatever the state.
sol.bprev = [];
for it=1:opt.maxit
    bprev = zeros(opt.nk, opt.ny);
    for i=1:opt.ny
        y = exp(sol.logy(i));
        % Next period log income at the quadrature nodes, for every b(t).
        ly1 = par.RHO*sol.logy(i)+par.SIGMA*e;
        B = repmat(bnext(:,i), 1, opt.nq);
        LY = repmat(ly1', opt.nk, 1);
        [~, c1] = buffer_policy(sol, B(:), LY(:));
        mu = reshape(c1.^(-par.GAMMAC), opt.nk, opt.nq)*w;
        c = (par.BETA*par.R*mu).^(-1/par.GAMMAC);
        bprev(:,i) = (y+bnext(:,i)-c)/par.R;
    end
    if ~isempty(sol.bprev)
        dist = max(abs(bprev(:)-sol.bprev(:)));
    else
        dist = Inf;
    end
    sol.bprev = bprev;
    if dist<opt.tol
        break
    end
end
sol.iterations = it;
sol.distance = dist;


function [x, w] = gauss_hermite(n)
% Nodes and weights for the expectation of a standard normal variable.
i = 1:n-1;
a = sqrt(i/2);
[V, D] = eig(diag(a, 1)+diag(a, -1));
[x, k] = sort(diag(D));
w = transpose(V(1,k).^2);
x = sqrt(2)*x;
