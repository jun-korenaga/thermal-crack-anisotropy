function [vsv,vsh,vpv,vph,xi,Bc,Bs,Gc,Gs,Hc,Hs,Ec,Es,A,L,F,N] = ...
    calc_vel_etc(C,rho)
% CALC_VEL_ETC Convert a 6-by-6 stiffness tensor in GPa to observables.

C = reshape(C,[6,6])*1e9; % now in Pa

vsv = sqrt(0.5*(C(4,4)+C(5,5))/rho);
vsh = sqrt(((C(1,1)+C(2,2))/8-C(1,2)/4+C(6,6)/2)/rho);
vpv = sqrt(C(3,3)/rho);
vph = sqrt(((C(1,1)+C(2,2))*(3/8)+C(1,2)/4+C(6,6)/2)/rho);
xi = (vsh/vsv)^2;

Bc = 0.5*(C(1,1)-C(2,2));
Bs = C(1,6)+C(2,6);
Gc = 0.5*(C(5,5)-C(4,4));
Gs = C(5,4);
Hc = 0.5*(C(1,3)-C(2,3));
Hs = C(3,6);
Ec = (C(1,1)+C(2,2))/8-C(1,2)/4-C(6,6)/2;
Es = (C(1,6)-C(2,6))/2;

A = 3/8*(C(1,1)+C(2,2))+1/4*C(1,2)+1/2*C(6,6);
L = 0.5*(C(4,4)+C(5,5));
F = 0.5*(C(1,3)+C(2,3));
N = 1/8*(C(1,1)+C(2,2))-1/4*C(1,2)+1/2*C(6,6);
end
