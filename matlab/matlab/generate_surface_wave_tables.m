function generate_surface_wave_tables()
% GENERATE_SURFACE_WAVE_TABLES Generate Rayleigh- and Love-wave data.

p = release_paths();
d0 = load_numeric_table(fullfile(p.prem,'Cij.iso.dat'),38);
r = d0(:,1); rho = d0(:,2); nr = length(r);
Ciso = tensors_from_table(d0);
istart = find(r>=(6371-410)*1e3,1);
w_full = trapezoid_weights(r);
w_upper = zeros(nr,1);
w_upper(istart:nr) = trapezoid_weights(r(istart:nr));

disp('precomputing Rayleigh-wave kernels');
R = rayleigh_kernels(p,r,rho,w_full,w_upper);
disp('precomputing Love-wave kernels');
L = love_kernels(p,r,rho,w_full,w_upper);

outdir = p.surface;
f_factor = 2;
R0 = R.R0_full;
L0 = L.L0_full;

write_ocean_references(p,R,L);

cases = {'1e-2',0.002;'1e-3',0.001};
maxlvls = [1 2 4 5];
single = {'','.single'};
lvl_rho = [2600 2900 3380.75 3377.73 3374.71];
dref0 = load_numeric_table(fullfile(p.prem,'Cij.prem.dat'),38);

for icase=1:size(cases,1)
  aspect = cases{icase,1}; phi_wanted = cases{icase,2};
  for im=1:length(maxlvls)
    maxlvl = maxlvls(im);
    for is=1:2
      dref = dref0;
      for ilvl=1:maxlvl
        name = ['Csat_prem' num2str(ilvl) '.a' aspect ...
                '.n3200' single{is} '.dat'];
        d = load_numeric_table(fullfile(p.tensors,name),38);
        [~,iphi] = min(abs(d(:,2)-phi_wanted));
        if abs(d(iphi,2)-phi_wanted)>10*eps(max(1,phi_wanted))
          error('porosity %.7g is absent from %s',phi_wanted,name);
        end
        Cij_mod = 1e9*d(iphi,3:38); % now in Pa
        range = dref(:,2)==lvl_rho(ilvl);
        if ~any(range)
          error('PREM layer density %.8g was not found',lvl_rho(ilvl));
        end
        dref(range,3:38) = repmat(Cij_mod,sum(range),1);
      end
      C = tensors_from_table(dref);
      coef = elastic_coefficients(C,Ciso);

      qR1 = coef.dA'*(R.wVA)+coef.dC'*(R.wVC) ...
          +f_factor*coef.dF'*(R.wVF)+coef.dL'*(R.wVL);
      qR2 = coef.Bc'*(R.wVA)+coef.Hc'*(2*R.wVF) ...
          +coef.Gc'*(R.wVL);
      qR4 = coef.Ec'*(R.wVA);
      dVR = (qR1'./R0)/(2)./(R.pvel*1e3)/1e3;
      dVRcos2 = (qR2'./R0)/(2)./(R.pvel*1e3)/1e3;
      dVRcos4 = (qR4'./R0)/(2)./(R.pvel*1e3)/1e3;
      tmp = [R.deg R.period R.pvel dVR dVRcos2 dVRcos4];
      tag = ['m' num2str(maxlvl) '.a' aspect '.p' ...
             num2str(phi_wanted) single{is} '.dat'];
      save(fullfile(outdir,['dVR.' tag]),'tmp','-ascii');

      qL1 = coef.dN'*(L.wWN)+coef.dL'*(L.wWL);
      qL2 = (-coef.Gc)'*(L.wWL);
      qL4 = (-coef.Ec)'*(L.wWN);
      dVL = (qL1'./L0)/(2)./(L.pvel*1e3)/1e3;
      dVLcos2 = (qL2'./L0)/(2)./(L.pvel*1e3)/1e3;
      dVLcos4 = (qL4'./L0)/(2)./(L.pvel*1e3)/1e3;
      tmp = [L.deg L.period L.pvel dVL dVLcos2 dVLcos4];
      save(fullfile(outdir,['dVL.' tag]),'tmp','-ascii');
    end
  end
end
end

function R = rayleigh_kernels(p,r,rho,w_full,w_upper)
d = load_numeric_table(fullfile(p.prem,'vel.prem_iso.3.dat'),4);
range = d(:,2)>=10 & d(:,2)<=250;
d = d(range,:); nL = size(d,1); nr = length(r);
R.deg = d(:,1); R.period = d(:,2); R.pvel = d(:,3);
R.wVA = zeros(nr,nL); R.wVC = zeros(nr,nL);
R.wVF = zeros(nr,nL); R.wVL = zeros(nr,nL);
R.R0_full = zeros(nL,1); R.R0_upper = zeros(nL,1);

for iL=1:nL
  L = R.deg(iL); k_sph = sqrt(L*(L+1));
  name = sprintf('S.0000000.%07d.ASC',L);
  d = load_numeric_table(fullfile(p.rayleigh,name),7);
  check_radius_grid(d(:,1),r,name);
  U = d(:,2); Udot = d(:,3); V = d(:,4); Vdot = d(:,5);
  [scale,deriv_scale] = eigenfunction_scales(R.period(iL));
  U = flipud(scale*U);
  dUdz = -flipud(deriv_scale*Udot);
  V = flipud(scale*k_sph*V);
  dVdz = -flipud(deriv_scale*k_sph*Vdot);
  k = 2*pi/(R.period(iL)*R.pvel(iL)*1e3);

  VA = V.^2; VC = (dUdz/k).^2;
  VF = dUdz.*V/k; VL = (dVdz/k-U).^2;
  R.wVA(:,iL) = w_upper.*VA;
  R.wVC(:,iL) = w_upper.*VC;
  R.wVF(:,iL) = w_upper.*VF;
  R.wVL(:,iL) = w_upper.*VL;
  energy = rho.*(U.^2+V.^2);
  R.R0_full(iL) = sum(w_full.*energy);
  R.R0_upper(iL) = sum(w_upper.*energy);
end
end

function Lout = love_kernels(p,r,rho,w_full,w_upper)
d = load_numeric_table(fullfile(p.prem,'vel.prem_iso.2.dat'),4);
range = d(:,2)>=10 & d(:,2)<=250;
d = d(range,:); nL = size(d,1); nr = length(r);
Lout.deg = d(:,1); Lout.period = d(:,2); Lout.pvel = d(:,3);
Lout.wWN = zeros(nr,nL); Lout.wWL = zeros(nr,nL);
Lout.L0_full = zeros(nL,1); Lout.L0_upper = zeros(nL,1);

for iL=1:nL
  L = Lout.deg(iL); k_sph = sqrt(L*(L+1));
  name = sprintf('T.0000000.%07d.ASC',L);
  d = load_numeric_table(fullfile(p.love,name),3);
  check_radius_grid(d(:,1),r,name);
  W = d(:,2); Wdot = d(:,3);
  W(1:2) = 0; Wdot(1:2) = 0; % zero out the ocean layer
  [scale,deriv_scale] = eigenfunction_scales(Lout.period(iL));
  W = flipud(scale*k_sph*W);
  dWdz = -flipud(deriv_scale*k_sph*Wdot);
  k = 2*pi/(Lout.period(iL)*Lout.pvel(iL)*1e3);

  WN = W.^2; WL = (dWdz/k).^2;
  Lout.wWN(:,iL) = w_upper.*WN;
  Lout.wWL(:,iL) = w_upper.*WL;
  energy = rho.*W.^2;
  Lout.L0_full(iL) = sum(w_full.*energy);
  Lout.L0_upper(iL) = sum(w_upper.*energy);
end
end

function [scale,deriv_scale] = eigenfunction_scales(period)
G = 6.6723e-11; Rmax = 6371e3; rho_scale = 5515;
vel_scale = Rmax*sqrt(pi*G*rho_scale);
omega_scale = vel_scale/Rmax;
omega_n = (2*pi/period)/omega_scale;
U_scale = sqrt(1/(rho_scale*Rmax^3));
scale = omega_n*U_scale;
deriv_scale = omega_n*U_scale/Rmax;
end

function check_radius_grid(r_eig,r,name)
if length(r_eig)~=length(r) || max(abs(flipud(r_eig)-r))>1e-6
  error('eigenfunction radius grid does not match Cij grid: %s',name);
end
end

function write_ocean_references(p,R,L)
d = load_numeric_table(fullfile(p.prem,'vel.prem_ocean.3.dat'),4);
[tf,loc] = ismember(R.deg,d(:,1));
if ~all(tf), error('Rayleigh ocean reference degrees do not match'); end
tmp = [d(loc,2) d(loc,3)];
save(fullfile(p.common,'ref.Rayleigh.dat'),'tmp','-ascii');

d = load_numeric_table(fullfile(p.prem,'vel.prem_ocean.2.dat'),4);
[tf,loc] = ismember(L.deg,d(:,1));
if ~all(tf), error('Love ocean reference degrees do not match'); end
tmp = [d(loc,2) d(loc,3)];
save(fullfile(p.common,'ref.Love.dat'),'tmp','-ascii');
end

function C = tensors_from_table(d)
nr = size(d,1); C = zeros(nr,6,6);
for i=1:nr
  C(i,:,:) = reshape(d(i,3:38),[6 6]);
end
end

function c = elastic_coefficients(C,C0)
nr = size(C,1);
names = {'dA','dC','dF','dL','dN','Bc','Hc','Gc','Ec'};
for iname=1:length(names)
  c.(names{iname}) = zeros(nr,1);
end
for i=1:nr
  Ct = reshape(C(i,:,:),[6 6]);
  C0t = reshape(C0(i,:,:),[6 6]); dC = Ct-C0t;
  c.dA(i) = 3/8*(dC(1,1)+dC(2,2))+1/4*dC(1,2)+1/2*dC(6,6);
  c.dC(i) = dC(3,3);
  c.dF(i) = 0.5*(dC(1,3)+dC(2,3));
  c.dL(i) = 0.5*(dC(4,4)+dC(5,5));
  c.dN(i) = 1/8*(dC(1,1)+dC(2,2))-1/4*dC(1,2)+1/2*dC(6,6);
  c.Bc(i) = 0.5*(Ct(1,1)-Ct(2,2));
  c.Hc(i) = 0.5*(Ct(1,3)-Ct(2,3));
  c.Gc(i) = 0.5*(Ct(5,5)-Ct(4,4));
  c.Ec(i) = 1/8*(Ct(1,1)+Ct(2,2))-1/4*Ct(1,2)-1/2*Ct(6,6);
end
end

function w = trapezoid_weights(x)
n = length(x); w = zeros(n,1);
if n<2, error('at least two radius samples are required'); end
w(1) = (x(2)-x(1))/2;
w(2:n-1) = (x(3:n)-x(1:n-2))/2;
w(n) = (x(n)-x(n-1))/2;
end
