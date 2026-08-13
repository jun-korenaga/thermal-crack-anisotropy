function generate_sks_tables()
% GENERATE_SKS_TABLES Generate maximum SKS delay curves for the main figure.

p = release_paths();
dref0 = load_numeric_table(fullfile(p.prem,'Cij.prem.dat'),38);
r = dref0(:,1); rho = dref0(:,2); nr = length(r);
istart = find(r>=(6371-410)*1e3,1);

aspect = {'1e-2','1e-3','1e-4'};
lvl_rho = [2600 2900 3380.75 3377.73 3374.71];
maxlvls = [1 2 4 5];

for im=1:length(maxlvls)
  maxlvl = maxlvls(im);
  for ia=1:3
    inputs = cell(maxlvl,1);
    for ilvl=1:maxlvl
      name = ['Csat_prem' num2str(ilvl) '.a' aspect{ia} ...
              '.n3200.single.dat'];
      inputs{ilvl} = load_numeric_table(fullfile(p.tensors,name),38);
    end
    nphi = size(inputs{1},1); out = zeros(nphi,2);
    for iphi=1:nphi
      dref = dref0;
      for ilvl=1:maxlvl
        Cij_mod = 1e9*inputs{ilvl}(iphi,3:38); % now in Pa
        range = dref(:,2)==lvl_rho(ilvl);
        dref(range,3:38) = repmat(Cij_mod,sum(range),1);
      end
      C = zeros(nr,6,6);
      for i=1:nr
        C(i,:,:) = reshape(dref(i,3:38),[6 6]);
      end
      vals = zeros(nr,1);
      for i=istart:nr
        vals(i) = int_SKS_deltat(i,rho,C);
      end
      delta_t = abs(trapz(r(istart:nr),vals(istart:nr)));
      out(iphi,:) = [100*inputs{1}(iphi,2) delta_t];
    end
    name = ['SKS_dtmax.a' num2str(ia) '.m' num2str(maxlvl) '.dat'];
    save(fullfile(p.common,name),'out','-ascii');
  end
end
end
