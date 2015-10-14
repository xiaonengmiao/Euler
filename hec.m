% Harten's entropy collection
% Author Luo Li, lluoac@ust.hk
function [f]=hec(x)
delta=200.;
f=max(abs(x),min(delta,0.5*(x*x+delta*delta)/delta));