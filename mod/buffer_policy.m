function [b, c] = buffer_policy(sol, bprev, logy)

% Evaluates the global solution of the borrowing constraint model.
%
% INPUTS
% - sol         [struct]    solution returned by buffer_egm
% - bprev       [double]    N×1 debt of the previous period
% - logy        [double]    N×1 log income of the current period
%
% OUTPUTS
% - b           [double]    N×1 debt of the current period
% - c           [double]    N×1 consumption

par = sol.par;
bprev = bprev(:);
logy = logy(:);
y = exp(logy);
limit = par.M*y;

if isempty(sol.bprev)
    b = limit;
else
    % Unconstrained policy at the two income nodes surrounding each point, linear
    % in log income in between (and linearly extrapolated outside the grid).
    ny = length(sol.logy);
    h = sol.logy(2)-sol.logy(1);
    j = min(max(floor((logy-sol.logy(1))/h)+1, 1), ny-1);
    t = (logy-sol.logy(j)')/h;
    bu = zeros(size(bprev));
    for i=unique(j)'
        k = j==i;
        lo = interp1(sol.bprev(:,i), sol.bnext(:,i), bprev(k), 'linear', 'extrap');
        hi = interp1(sol.bprev(:,i+1), sol.bnext(:,i+1), bprev(k), 'linear', 'extrap');
        bu(k) = (1-t(k)).*lo+t(k).*hi;
    end
    b = min(bu, limit);
end

c = y+b-par.R*bprev;
