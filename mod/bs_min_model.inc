// The borrowing constraint model with the complementarity written as a single
// equation, lb >= 0, M*y-b >= 0 and one of them zero: the formulation the
// extended path can solve, and the yardstick of the comparison.

@#include "bs_declarations.inc"

model;
  [name='expected consumption']
  ec = c(1);
  [name='budget constraint']
  c = y + b - R*b(-1);
  [name='borrowing']
  min(lb, M*y - b) = 0;
  [name='euler']
  lb = 1/c^GAMMAC - BETA*R/c(+1)^GAMMAC;
  [name='income']
  log(y) = RHO*log(y(-1)) + u;
end;

@#include "bs_calibration.inc"

options_.nograph = true;
options_.ep.periods = 100;
