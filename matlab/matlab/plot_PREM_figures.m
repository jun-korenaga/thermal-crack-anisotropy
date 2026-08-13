% plot_PREM_figures.m
% make the PREM1 and PREM2 surface-wave figures

function plot_PREM_figures(datadir,figdir)

if nargin<1
  p = release_paths();
  datadir = p.results;
end
if nargin<2
  p = release_paths();
  figdir = p.figures;
end

make_plot('PREM1','a1e-2.p0.002',datadir,figdir);
make_plot('PREM2','a1e-3.p0.001',datadir,figdir);
end

function make_plot(id,id_data,datadir,figdir)

dir = fullfile(datadir,'surface_waves');
dir0 = fullfile(datadir,'common');
refR = load(fullfile(dir0,'ref.Rayleigh.dat'));
refL = load(fullfile(dir0,'ref.Love.dat'));

maxlvl = [1 2 4 5];
linec = [141 160 203; 102 194 165; 231 138 195; 191 129 103]/255;

% axis limits and major ticks
if strcmp(id,'PREM1')
  ylim_all = {[-2.5 0],[0 1.6],[0 0.17],[-5 0],[-0.45 0],[-2 0]};
  yticks_all = {-2.5:0.5:0,0:0.5:1.5,0:0.05:0.15,...
                -5:1:0,-0.4:0.1:0,-2:0.5:0};
else
  ylim_all = {[-7 0],[0 5],[0 0.8],[-18 0],[-1.6 0],[-8 0]};
  yticks_all = {-7:1:0,0:1:5,0:0.2:0.8,...
                -15:5:0,-1.5:0.5:0,-8:2:0};
end

fig = figure('Visible','off','Color','w','Units','inches',...
             'Position',[1 1 7.79 5.33]);
set(fig,'PaperUnits','inches','PaperPosition',[0 0 7.79 5.33],...
        'PaperSize',[7.79 5.33]);

left = [0.08 0.405 0.73];
bottom = [0.56 0.10];
ax = gobjects(6,1);
for irow=1:2
  for icol=1:3
    ipanel = (irow-1)*3+icol;
    ax(ipanel) = axes(fig,'Position',[left(icol) bottom(irow) 0.26 0.37]);
    hold(ax(ipanel),'on'); box(ax(ipanel),'on');
    set(ax(ipanel),'FontName','Helvetica','FontSize',9,...
        'LineWidth',0.9,'TickDir','in','XMinorTick','on','YMinorTick','on',...
        'XLim',[10 250],'XTick',50:50:250,...
        'YLim',ylim_all{ipanel},'YTick',yticks_all{ipanel},...
        'Layer','top');
  end
end

% plot the four assumed crack depths
for im=4:-1:1
  level = maxlvl(im);
  tag = ['m' num2str(level) '.' id_data];
  dR_single = load(fullfile(dir,['dVR.' tag '.single.dat']));
  dR_hex = load(fullfile(dir,['dVR.' tag '.dat']));
  dL_single = load(fullfile(dir,['dVL.' tag '.single.dat']));
  dL_hex = load(fullfile(dir,['dVL.' tag '.dat']));

  check_periods(refR,dR_single,['Rayleigh ' tag]);
  check_periods(refL,dL_single,['Love ' tag]);

  vani = refR(:,2); viso = dR_single(:,3);
  val = (dR_single(:,4)-(vani-viso))./vani*100;
  plot(ax(1),refR(:,1),val,'--','Color',linec(im,:),'LineWidth',1.3);
  viso = dR_hex(:,3);
  val = (dR_hex(:,4)-(vani-viso))./vani*100;
  plot(ax(1),refR(:,1),val,'-','Color',linec(im,:),'LineWidth',1.3);
  plot(ax(2),refR(:,1),dR_single(:,5)./vani*100,'--',...
       'Color',linec(im,:),'LineWidth',1.3);
  plot(ax(3),refR(:,1),dR_single(:,6)./vani*100,'--',...
       'Color',linec(im,:),'LineWidth',1.3);

  vani = refL(:,2); viso = dL_single(:,3);
  val = (dL_single(:,4)-(vani-viso))./vani*100;
  plot(ax(4),refL(:,1),val,'--','Color',linec(im,:),'LineWidth',1.3);
  viso = dL_hex(:,3);
  val = (dL_hex(:,4)-(vani-viso))./vani*100;
  plot(ax(4),refL(:,1),val,'-','Color',linec(im,:),'LineWidth',1.3);
  plot(ax(5),refL(:,1),dL_single(:,5)./vani*100,'--',...
       'Color',linec(im,:),'LineWidth',1.3);
  plot(ax(6),refL(:,1),dL_single(:,6)./vani*100,'--',...
       'Color',linec(im,:),'LineWidth',1.3);
end

title(ax(1),'constant','FontWeight','normal','FontSize',12);
title(ax(2),'cos 2\psi','FontWeight','normal','FontSize',12);
title(ax(3),'cos 4\psi','FontWeight','normal','FontSize',12);
ylabel(ax(1),'Rayleigh \delta c/c [%]','FontSize',11);
ylabel(ax(4),'Love \delta c/c [%]','FontSize',11);
for i=4:6
  xlabel(ax(i),'period [s]','FontSize',11);
end

% panel labels
panel = {'(a)','(b)','(c)','(d)','(e)','(f)'};
text(ax(1),0.92,0.08,panel{1},'Units','normalized',...
     'HorizontalAlignment','right','FontSize',11);
for i=2:3
  text(ax(i),0.92,0.90,panel{i},'Units','normalized',...
       'HorizontalAlignment','right','FontSize',11);
end
for i=4:6
  text(ax(i),0.92,0.08,panel{i},'Units','normalized',...
       'HorizontalAlignment','right','FontSize',11);
end

% legend in panel (e)
h = gobjects(6,1);
for im=1:4
  h(im) = plot(ax(5),nan,nan,'-','Color',linec(im,:),'LineWidth',1.3);
end
h(5) = plot(ax(5),nan,nan,'k-','LineWidth',1.3);
h(6) = plot(ax(5),nan,nan,'k--','LineWidth',1.3);
leg = legend(ax(5),h,{'d=15 km','d=25 km','d=52 km','d=80 km',...
             'hexagonal','single'},'Location','south','FontSize',8.5);
set(leg,'Box','on','Color','w','EdgeColor','k');

drawnow;
epsfile = fullfile(figdir,[id '.eps']);
print(fig,epsfile,'-depsc2','-vector','-r600');
copyfile(epsfile,fullfile(figdir,[id '.ps']),'f');
close(fig);
end

function check_periods(ref,d,name)
if size(ref,1)~=size(d,1)
  error('row-count mismatch for %s',name);
end
end
