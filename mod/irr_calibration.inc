// Calibration shared by the two formulations. Taken from
// tests/occbin/model_irrcap_twoconstraints/dynrbc_common.inc.

steady_state_model;
  kss = ((1/BETA-1+DELTAK)/ALPHA)^(1/(ALPHA-1));
  css = -DELTAK*kss + kss^ALPHA;
  iss = DELTAK*kss;
  k = log(kss);
  c = log(css);
  i = log(iss);
  lambdak = 0;
  a = 0;
end;

BETA   = 0.96;
ALPHA  = 0.33;
DELTAK = 0.10;
GAMMAC = 2;
RHOA   = 0.9;
PHI    = 0.975;

shocks;
  var erra; stderr 0.015;
end;

steady;
