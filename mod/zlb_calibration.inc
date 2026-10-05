// Calibration of tests/ep/hybrid/rs.mod (Rudebusch and Swanson, 2012), but for
// the preferences and the discount factor. The risk aversion (CRRA 20) and the
// intertemporal elasticity of substitution (0.5) are milder than in the paper
// (110 and 0.11): with the paper's values the certainty equivalent raises V to a
// power close to 50, and the extended path solvers fail away from the steady
// state. The nominal rate is 2% a year at the steady state, so that the zero
// lower bound binds in a sizeable fraction of the periods.
//
// The macro variables below can be set before the inclusion; a shock with a zero
// standard error is neither simulated nor integrated over.

@#ifndef IES
  @#define IES = 0.5
@#endif
@#ifndef CRRA
  @#define CRRA = 20
@#endif
@#ifndef DZBAR
  @#define DZBAR = 1.0025
@#endif
@#ifndef RBETA
  @#define RBETA = 0.995
@#endif
@#ifndef RHOA
  @#define RHOA = 0.96
@#endif
@#ifndef RHOZ
  @#define RHOZ = 0.85
@#endif
@#ifndef SIGA
  @#define SIGA = 0.002
@#endif
@#ifndef SIGZ
  @#define SIGZ = 0.005
@#endif
@#ifndef SIGINT
  @#define SIGINT = 0
@#endif

DZBar = @{DZBAR};
eta = 2/3;
phi = 1/@{IES};
beta = @{RBETA}*DZBar^phi;
LMax = 3;
chi0 = 1/3;
chi = 1/.28*(LMax-1);
alpha = (@{CRRA} - 1/(1/phi + 1/chi*(LMax-1)))*(1/(1-phi) + 1/(1-chi)*(LMax-1));
theta = .2;
xi = .78;
YBar = 10^((1-eta)/eta);
GBar = .17*YBar;
KBar = 10*YBar;
IBar = (.02+(DZBar-1))*KBar;
piBar = 1.0;
taylrho = .73;
taylpi = .53;
tayly = .93;
rhoa = @{RHOA};
rhoz = @{RHOZ};
rhoinflavg = .7;
VAIMSS = 1;

steady_state_model;
  lA = 0;
  lz = 0;
  lpi = log(piBar);
  lpiavg = lpi;
  lDisp = 0;
  Int = log(DZBar^phi/beta)+lpi;
  Ishadow = Int;
  lL = 0;
  lY = log(YBar);
  lC = log(exp(lY)-GBar-IBar);
  lp0 = log(((exp(lpi)^(-1/theta)-xi)/((1-xi)*exp(lpi)^(-1/theta)))^(-theta));
  lzd = log(exp(lY)/(1-xi*beta*DZBar^(-phi)*exp(lpi)^(1/theta)));
  lzn = log(exp(lp0)^(1+(1+theta)/theta*(1-eta)/eta)*exp(lzd));
  MC = (exp(lzn)-xi*beta*DZBar^(-phi)*piBar^((1+theta)/theta/eta)*exp(lzn))/((1+theta)*YBar);
  lwreal = log(MC*eta/(exp(lY)^((1-eta)/eta)/KBar^((1-eta)/eta)));
  chi0 = exp(lwreal)*exp(lC)^(-phi)/(LMax-exp(lL))^(-chi);
  V = (exp(lC)^(1-phi)/(1-phi) + chi0*(LMax-exp(lL))^(1-chi)/(1-chi))/(1-beta*DZBar^(1-phi));
  VAIMSS = 1.05*V;
  Valphaexp = (V/VAIMSS)^(1-alpha);
  Vkp = DZBar^(1-phi)*V;
end;

shocks;
  var epsA; stderr @{SIGA};
  var epsZ; stderr @{SIGZ};
  var epsInt; stderr @{SIGINT};
end;

steady;
