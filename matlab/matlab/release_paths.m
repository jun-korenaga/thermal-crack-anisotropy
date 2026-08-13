function p = release_paths()
% RELEASE_PATHS Return paths relative to this package.

p.matlab = fileparts(mfilename('fullpath'));
p.root = fileparts(p.matlab);
p.data = fullfile(p.root,'data');
p.tensors = fullfile(p.data,'effective_tensors');
p.prem = fullfile(p.data,'prem');
p.rayleigh = fullfile(p.prem,'eigen_rayleigh');
p.love = fullfile(p.prem,'eigen_love');
p.results = fullfile(p.root,'results');
p.common = fullfile(p.results,'common');
p.surface = fullfile(p.results,'surface_waves');
p.figures = fullfile(p.root,'figures');
end
