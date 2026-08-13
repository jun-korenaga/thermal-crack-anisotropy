function generate_velocity_tables()
% GENERATE_VELOCITY_TABLES Generate inputs for manuscript CPO/TC figures.

p = release_paths();
id = {'jungA','abram97'};
rho = [3320 3355];
aspect = {'1e-2','1e-3','1e-4'};
single = {'','.single'};

for iid=1:length(id)
  for is=1:2
    for ia=1:3
      name = ['Csat_' id{iid} '.a' aspect{ia} '.n3200' single{is} '.dat'];
      process_one(fullfile(p.tensors,name),rho(iid),...
                  fullfile(p.common,['vel_' id{iid} '.a' aspect{ia} ...
                  '.n3200' single{is} '.dat']));
    end
  end
end

% Longer single-crack sequences used in the scaling figure.
for ia=1:3
  name = ['Csat_jungA.a' aspect{ia} '.n3200.single.dat1c'];
  process_one(fullfile(p.tensors,name),3320,...
              fullfile(p.common,['vel_jungA.a' aspect{ia} ...
              '.n3200.single.dat1c']));
end
end

function process_one(infile,rho,outfile)
d = load_numeric_table(infile,38);
phi = d(:,2); Cij = d(:,3:38); n = length(phi);
tmp = zeros(n,18); tmp(:,1) = phi;

for i=1:n
  [vsv,vsh,vpv,vph,xi,Bc,Bs,Gc,Gs,Hc,Hs,Ec,Es,A,L,F,N] = ...
      calc_vel_etc(Cij(i,:),rho);
  tmp(i,:) = [phi(i) vsv vsh xi vpv vph Bc Bs Gc Gs Hc Hs ...
              Ec Es A L F N];
end
save(outfile,'tmp','-ascii');
end
