% Roe flux for euler equation
% ul - 1x3 array, state on the left
% ur - 1x3 array, state on the right
% flux - 1x1x3 array, output flux
% Author Luo Li, lluoac@ust.hk
function [flux]=r_f_euler(ul,ur)
gammer=1.4;
flux=zeros(1,3);

% left state
rhol=ul(1);
rhoul=ul(2);
rhoEl=ul(3);
% right state
rhor=ur(1);
rhour=ur(2);
rhoEr=ur(3);
	    
% primitive variables
ulx=rhoul/rhol;
urx=rhour/rhor;
pl=(gammer-1.)*(rhoEl-0.5*rhoul*rhoul/rhol);
pr=(gammer-1.)*(rhoEr-0.5*rhour*rhour/rhor);
hl=(rhoEl+pl)/rhol;
hr=(rhoEr+pr)/rhor;

% Roe averages
rr=sqrt(rhol/rhor);
rhoavg=rr*rhor;
uavg=(rr*ulx+urx)/(rr+1.);
havg=(rr*hl+hr)/(rr+1.);
aavg=sqrt((gammer-1.)*(havg-0.5*uavg*uavg));

% wavespeeds
lambda0=uavg;
lambdap=uavg+aavg;
lambdam=uavg-aavg;

% wave strenghs
dv0=rhol-rhor-(pl-pr)/(aavg*aavg);
dvp=ulx-urx+(pl-pr)/(rhoavg*aavg);
dvm=ulx-urx-(pl-pr)/(rhoavg*aavg);

% flux
flux(1,1)=0.5*(rhoul+rhour)-0.5*hec(lambda0)*dv0-...
			0.25*hec(lambdap)*dvp*rhoavg/aavg+...
			0.25*hec(lambdam)*dvm*rhoavg/aavg;
flux(1,2)=0.5*(rhoul*ulx+pl+rhour*urx+pr)-0.5*hec(lambda0)*dv0*uavg-...
			0.25*hec(lambdap)*dvp*rhoavg*lambdap/aavg+...
			0.25*hec(lambdam)*dvm*rhoavg*lambdam/aavg;
flux(1,3)=0.5*(rhoEl*ulx+pl*ulx+rhoEr*urx+pr*urx)-...
			0.25*hec(lambda0)*dv0*uavg*uavg-...
			0.25*hec(lambdap)*dvp*rhoavg*(havg+aavg*uavg)/aavg+...
			0.25*hec(lambdam)*dvm*rhoavg*(havg-aavg*uavg)/aavg;

