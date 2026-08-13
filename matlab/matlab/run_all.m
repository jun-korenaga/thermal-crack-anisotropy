function run_all()
% RUN_ALL Generate all numerical tables used for the main figures.

p = release_paths();
make_output_dirs(p);
self_test();

disp('generating effective-velocity tables');
generate_velocity_tables();
disp('generating scaling curves');
generate_tcscale3_tables();
disp('generating SKS curves');
generate_sks_tables();
disp('generating surface-wave tables');
generate_surface_wave_tables();

validate_release();
disp(['finished: results are under ' p.results]);
end

function make_output_dirs(p)
dirs = {p.results,p.common,p.surface,p.figures};
for i=1:length(dirs)
  if ~exist(dirs{i},'dir')
    mkdir(dirs{i});
  end
end
end
