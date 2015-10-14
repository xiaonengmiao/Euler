% Kinetic flux for euler equation
% wl - states on the left
% wr - states on the right
% flux - output flux
% Author Haohan Li, hlibb@connect.ust.hk
function [flux] = k_f_Euler(wl, wr)
flux = zeros(1, 3);
gammar = 1.4;
K = 4;
pi = 3.1415926;

%% left states
rhol = wl(1);
rhoul = wl(2);
rhoel = wl(3);
ul = rhoul/rhol;
pl = (gammar-1.0)*(rhoel-0.5*rhol*ul*ul);
laml = rhol/(2.0*pl);
erl = erfc(-sqrt(laml)*ul);

%% left states
rhor = wr(1);
rhour = wr(2);
rhoer = wr(3);
ur = rhour/rhor;
pr = (gammar-1.0)*(rhoer-0.5*rhor*ur*ur);
lamr = rhor/(2.0*pr);
err = erfc(sqrt(lamr)*ur);

%% flux
flux(1) = rhol*(0.5*ul*erl+0.5*exp(-laml*ul*ul)/sqrt(pi*laml))*1....
    +rhor*(0.5*ur*err-0.5*exp(-lamr*ur*ur)/sqrt(pi*lamr))*1.;
flux(2) = rhol*((0.5*ul*ul+1/(4*laml))*erl+0.5*ul*exp(-laml*ul*ul)/sqrt(pi*laml))...
    +rhor*((0.5*ur*ur+1/(4*lamr))*err-0.5*ur*exp(-lamr*ur*ur)/sqrt(pi*lamr));
flux(3) = rhol*((1/4*ul*ul*ul+(K+3)/(8*laml)*ul)*erl+(1/4*ul*ul+(K+2)/(8*laml))*exp(-laml*ul*ul)/sqrt(pi*laml))...
    +rhor*((1/4*ur*ur*ur+(K+3)/(8*lamr)*ur)*err-(1/4*ur*ur+(K+2)/(8*lamr))*exp(-lamr*ur*ur)/sqrt(pi*lamr));