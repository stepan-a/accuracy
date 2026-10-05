// Equations of the New Keynesian model with Epstein-Zin preferences, but for the
// zero lower bound, which is written in zlb_max_model.inc, zlb_mcp_model.inc and
// occbin_zlb.mod. Compared with tests/ep/hybrid/rs.mod: no trend growth shock nor
// inflation target shock, no financial variables, a discount factor shock lz,
// government spending held at its steady state, a shadow rate Ishadow set by the
// Taylor rule, and none of the bounds on Valphaexp and lL.

  [name='value']
  V = exp(lC)^(1-phi)/(1-phi) + chi0*(LMax-exp(lL))^(1-chi)/(1-chi) + beta*exp(lz)*Vkp;
  [name='certainty equivalent']
  Valphaexp = (V(+1)/VAIMSS)^(1-alpha);
  [name='risk adjusted value']
  Vkp = VAIMSS*DZBar^(1-phi)*Valphaexp^(1/(1-alpha));
  [name='euler']
  exp(lC)^(-phi) = beta*exp(lz)*exp(Int-lpi(+1))*exp(lC(+1))^(-phi)*DZBar^(-phi)
                   *(V(+1)*DZBar^(1-phi)/Vkp)^(-alpha);
  [name='price numerator']
  exp(lzn) = (1+theta)*MC*exp(lY) + xi*beta*exp(lz)*exp(lC(+1)-lC)^(-phi)*DZBar^(-phi)
             *(V(+1)*DZBar^(1-phi)/Vkp)^(-alpha)*exp(lpi(+1))^((1+theta)/theta/eta)*exp(lzn(+1));
  [name='price denominator']
  exp(lzd) = exp(lY) + xi*beta*exp(lz)*exp(lC(+1)-lC)^(-phi)*DZBar^(-phi)
             *(V(+1)*DZBar^(1-phi)/Vkp)^(-alpha)*exp(lpi(+1))^(1/theta)*exp(lzd(+1));
  [name='reset price']
  exp(lp0)^(1+(1+theta)/theta*(1-eta)/eta) = exp(lzn-lzd);
  [name='inflation']
  exp(lpi)^(-1/theta) = (1-xi)*exp(lp0+lpi)^(-1/theta) + xi;
  [name='marginal cost']
  MC = exp(lwreal)/eta*exp(lY)^((1-eta)/eta)/exp(lA)^(1/eta)/KBar^((1-eta)/eta);
  [name='labour supply']
  chi0*(LMax-exp(lL))^(-chi)/exp(lC)^(-phi) = exp(lwreal);
  [name='production']
  exp(lY) = exp(lA)*KBar^(1-eta)*exp(lL)^eta/exp(lDisp);
  [name='price dispersion']
  exp(lDisp)^(1/eta) = (1-xi)*exp(lp0)^(-(1+theta)/theta/eta) + xi*exp(lpi)^((1+theta)/theta/eta)*exp(lDisp(-1))^(1/eta);
  [name='resource constraint']
  exp(lC) = exp(lY) - GBar - IBar;
  [name='average inflation']
  lpiavg = rhoinflavg*lpiavg(-1) + (1-rhoinflavg)*lpi;
  [name='taylor rule']
  4*Ishadow = (1-taylrho)*(4*log(DZBar^phi/beta) + 4*lpiavg + taylpi*(4*lpiavg-4*log(piBar)) + tayly*(exp(lY)-YBar)/YBar)
              + taylrho*4*Int(-1) + epsInt;
  [name='technology']
  lA = rhoa*lA(-1) + epsA;
  [name='discount factor']
  lz = rhoz*lz(-1) + epsZ;
