// Declarations shared by the OccBin and the extended path formulations of the
// New Keynesian model with Epstein-Zin preferences and a zero lower bound on the
// nominal interest rate, so that they are the same economy. Streamlined from
// tests/ep/hybrid/rs.mod (Rudebusch and Swanson, 2012).

var V Vkp Valphaexp lC lL lY lpi lp0 lzn lzd MC lwreal lDisp lpiavg Int Ishadow lA lz;

varexo epsA epsZ epsInt;

parameters theta xi beta phi alpha eta KBar chi0 LMax chi IBar GBar YBar DZBar piBar
           rhoinflavg taylrho taylpi tayly rhoa rhoz VAIMSS;
