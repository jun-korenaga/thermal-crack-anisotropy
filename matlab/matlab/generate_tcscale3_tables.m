function generate_tcscale3_tables()
% GENERATE_TCSCALE3_TABLES Generate curves used in manuscript Figure TCscale3.

p = release_paths();
aspect = {'1e-2','1e-3','1e-4'};
d = cell(3,1);
for ia=1:3
  name = ['vel_jungA.a' aspect{ia} '.n3200.single.dat1c'];
  d{ia} = load_numeric_table(fullfile(p.common,name),18);
end

rho = 3320; is = 2;
xi0 = 1.07; B_A0 = 0.03; G_L0 = 0.02; E_N0 = 0.01;
contour_BA = []; contour_GL = []; contour_EN = []; contour_xi = [];

for ia=1:3
  [new_xi,new_B_A,new_G_L,new_E_N] = ...
      apply_scaling(is,d{ia},rho,xi0,B_A0,G_L0,E_N0);
  xvec = {new_B_A*100,new_G_L*100,new_E_N*100};
  for i=1:3
    tmp = [xvec{i} new_xi];
    save(fullfile(p.common,['TCscale3.a' num2str(ia) '.i' ...
         num2str(i) '.dat']),'tmp','-ascii');
  end
  contour_BA = [contour_BA; xvec{1}']; %#ok<AGROW>
  contour_GL = [contour_GL; xvec{2}']; %#ok<AGROW>
  contour_EN = [contour_EN; xvec{3}']; %#ok<AGROW>
  contour_xi = [contour_xi; new_xi']; %#ok<AGROW>
end

save(fullfile(p.common,'TCscale3.con_BA.dat'),'contour_BA','-ascii');
save(fullfile(p.common,'TCscale3.con_GL.dat'),'contour_GL','-ascii');
save(fullfile(p.common,'TCscale3.con_EN.dat'),'contour_EN','-ascii');
save(fullfile(p.common,'TCscale3.con_xi.dat'),'contour_xi','-ascii');
end
