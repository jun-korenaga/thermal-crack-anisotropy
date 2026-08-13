function validate_release()
% VALIDATE_RELEASE Check the expected outputs and their basic invariants.

p = release_paths();
aspect = {'1e-2','1e-3','1e-4'};
single = {'','.single'};

for iid={'jungA','abram97'}
  for ia=1:3
    for is=1:2
      name = ['vel_' iid{1} '.a' aspect{ia} '.n3200' single{is} '.dat'];
      d = load_numeric_table(fullfile(p.common,name),18);
      assert(all(diff(d(:,1))>=0),'porosity must be monotonic');
      assert(all(d(:,2:6)>0,'all'),'wave speeds must be positive');
    end
  end
end

for maxlvl=[1 2 4 5]
  for ia=1:3
    d = load_numeric_table(fullfile(p.common,['SKS_dtmax.a' ...
        num2str(ia) '.m' num2str(maxlvl) '.dat']),2);
    assert(all(d(:,2)>=0),'SKS delays must be nonnegative');
  end
end

cases = {'1e-2','0.002';'1e-3','0.001'};
for icase=1:size(cases,1)
  for maxlvl=[1 2 4 5]
    for is=1:2
      suffix = ['m' num2str(maxlvl) '.a' cases{icase,1} ...
                '.p' cases{icase,2} single{is} '.dat'];
      dR = load_numeric_table(fullfile(p.surface,['dVR.' suffix]),6);
      dL = load_numeric_table(fullfile(p.surface,['dVL.' suffix]),6);
      assert(all(dR(:,2)>=10 & dR(:,2)<=250),'Rayleigh period out of range');
      assert(all(dL(:,2)>=10 & dL(:,2)<=250),'Love period out of range');
    end
  end
end
disp('validation passed');
end
