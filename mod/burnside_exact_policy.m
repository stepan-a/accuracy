function Y = burnside_exact_policy(S, U, M_, oo_)

% Exact solution of the model of Burnside (1998).
%
% INPUTS
% - S           [double]    1×N lagged value of x, in deviation from the
%                           deterministic steady state
% - U           [double]    1×N innovations
% - M_          [struct]    describing the model
% - oo_         [struct]    results
%
% OUTPUTS
% - Y           [double]    2×N current endogenous variables, in declaration order
%
% REMARKS
% y(t) = sum_{i>=1} beta^i exp(a_i + b_i (x(t)-xbar)), with the coefficients
% below. They are copied from tests/decision_rules/second_order/burnside_1.mod,
% which checks them against the decision rule computed by Dynare.

% Copyright © 2026 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <https://www.gnu.org/licenses/>.

beta = M_.params(strcmp(M_.param_names, 'beta'));
theta = M_.params(strcmp(M_.param_names, 'theta'));
rho = M_.params(strcmp(M_.param_names, 'rho'));
xbar = M_.params(strcmp(M_.param_names, 'xbar'));
sigma2 = M_.Sigma_e(1,1);

i = 1:800;
a = theta*xbar*i ...
    + (theta^2*sigma2)/(2*(1-rho)^2) ...
      *(i - 2*rho*(1-rho.^i)/(1-rho) + rho^2*(1-rho.^(2*i))/(1-rho^2));
b = theta*rho*(1-rho.^i)/(1-rho);
c = beta.^i.*exp(a);

ix = strcmp(M_.endo_names, 'x');
iy = strcmp(M_.endo_names, 'y');

% S is the lagged state in deviation from the deterministic steady state of x,
% which is xbar.
xlag = oo_.dr.ys(ix)+S;
x = (1-rho)*xbar+rho*xlag+U;

y = c*exp(transpose(b)*(x-xbar));

Y = zeros(M_.endo_nbr, size(S,2));
Y(ix,:) = x;
Y(iy,:) = y;
