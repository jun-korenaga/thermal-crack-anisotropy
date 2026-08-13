function [new_xi,new_B_A,new_G_L,new_E_N] = ...
    apply_scaling(is,data,rho,xi0,B_A0,G_L0,E_N0)
% APPLY_SCALING Apply thermal-crack changes to prescribed background CPO.
% is = 1 for hexagonal cracks and 2 for a single vertical crack set.

ivsv = 2; ivsh = 3; ixi = 4; ivph = 6;
iBc = 7; iGc = 9; iEc = 13;

xi_ref = data(1,ixi);
dxi_norm = (data(:,ixi)-xi_ref)/xi_ref;
new_xi = xi0*(1+dxi_norm);

A = rho*data(:,ivph).^2;
L = rho*data(:,ivsv).^2;
N = rho*data(:,ivsh).^2;
B = data(:,iBc); G = data(:,iGc); E = data(:,iEc);
dA_norm = (A-A(1))/A(1);
dL_norm = (L-L(1))/L(1);
dN_norm = (N-N(1))/N(1);

if is==1
  new_B_A = B_A0*(B/B(1))./(1+dA_norm);
  new_G_L = G_L0*(G/G(1))./(1+dL_norm);
  new_E_N = E_N0*(E/E(1))./(1+dN_norm);
else
  new_B_A = B_A0./(1+dA_norm)+(B-B(1))./A;
  new_G_L = G_L0./(1+dL_norm)+(G-G(1))./L;
  new_E_N = E_N0./(1+dN_norm)+(E-E(1))./N;
end
end
