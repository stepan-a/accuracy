// Calibration of tests/occbin/model_borrcon. With BETA*R<1 the household is
// impatient, and the borrowing limit binds at the deterministic steady state.

steady_state_model;
  b = M;
  c = 1+M-R*M;
  ec = c;
  lb = (1-BETA*R)/c^GAMMAC;
  y = 1;
end;

R = 1.05;
BETA = 0.945;
RHO = 0.9;
SIGMA = 0.05;
M = 1;
GAMMAC = 1;

shocks;
  var u; stderr SIGMA;
end;

steady;
