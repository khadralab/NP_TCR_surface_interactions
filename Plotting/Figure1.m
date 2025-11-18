%% Manuscript Figure 2
% 1) Comparison of different surfaces
addpath('Plotting/Functions')
addpath('Functions/')

%%
clear all;

f = figure('Color','white','units','centimeters','Position',[55 1 20 24],...
    'PaperUnits','centimeters','PaperPosition', [1 1 20 24], 'Renderer','painters');

half_hor = [0.52 0.67 0.82];
third_hor = [0.05, 0.2, 0.35];
third_ver = [0.9, 0.8, 0.7];
bot = [0.55, 0.45, 0.35];
left = [0.5, 0.65, 0.8];

for i=1:3
    ax2(i) = axes(f, 'Units','normalized','Position',[left(i)+0.05 0.88 0.1 0.08]); % TCR Surfaces
    ax5(i) = axes(f, 'Units','normalized','Position',[0.58 bot(i) 0.35 0.06]); % TCR coverage

    if i==3
        ax1 = axes(f, 'Units','normalized','Position',[0.1 third_ver(3) 0.2 0.28]); % NP Properties
        ax3 = axes(f, 'Units','normalized','Position',[0.6 third_ver(3) 0.35 0.28]); % Schematic Contact Area
        ax4 = axes(f, 'Units','normalized','Position',[0.1 0.35 0.35 0.22]); % NP capacity

        ax6 = axes(f, 'Units','normalized','Position',[0.1 0.05 0.35 0.22]); % Time series
        ax7 = axes(f, 'Units','normalized','Position',[0.6 0.05 0.35 0.22]); % Bound TCR Distributions
        break
    end

end

% New Simulations
koff_vals = [0.0001, 0.001, 0.01];
kd_vals = koff_vals ./ 0.1;
rho_vals = linspace(-4,3.5,16);
bound_type = {'phos', 'bound', 'np', 'covered','fracphos', 'normBound', 'normPhos','rescaled'};

valence = [5];
tpc = [1, 3, 20];
kp = [0.01,0];
np_radius = 20;

% Color and Linestyle
col = gray(length(valence)+2);
hol = winter(length(valence)+2);
sol = autumn(length(valence)+2);
pol = pink(length(valence)+2);

linestyle = {'-','--'};
markerstyle = {'s','d'};
colorstyle = {col, hol, sol, pol};

axis([ax1 ax3],'off');

% --------------------------------------------------------------------------
% NP Properties Schematic
%--------------------------------------------------------------------------

imshow('NP Manuscript Schematics.png','Parent',ax1)
set(ax1, 'Position',[0.02 0.56 0.45 0.48]);

% --------------------------------------------------------------------------
% TCR Surfaces
%--------------------------------------------------------------------------

a2_pos = get(ax2, 'Position');

rSurf = 1000;
tcr_per_cluster = 1;
num_tcr = 300;
max_clusters = floor(num_tcr/tcr_per_cluster);
num_clusters = [5,10,20];
cluster_radius = 50;
rTCR = 5;
rNP = 20;

col = [[0, 0, 0]; [0,0,1]; [1, 0, 0]];

tcr_pos = gen_clusters(rSurf, cluster_radius, 20, 0, num_tcr);
PlotSurface(ax2(1), rSurf, tcr_pos, rTCR, [], rNP, col(1,:));

tcr_pos = gen_clusters(rSurf, cluster_radius/sqrt(6), 100, 3, num_tcr);
PlotSurface(ax2(2), rSurf, tcr_pos, rTCR, [], rNP, col(2,:));

tcr_pos = gen_clusters(rSurf, cluster_radius, 15, 20, num_tcr);
PlotSurface(ax2(3), rSurf, tcr_pos, rTCR, [], rNP, col(3,:));

text(ax2(1),0.1, 2.2, '1 TPC','Units','centimeters','FontSize',14,'FontWeight','bold','Color',col(1,:));
text(ax2(2),0.1, 2.2, '3 TPC','Units','centimeters','FontSize',14,'FontWeight','bold','Color',col(2,:));
text(ax2(3),0, 2.2, '20 TPC','Units','centimeters','FontSize',14,'FontWeight','bold','Color',col(3,:));

set(ax2(1), 'Position', a2_pos{1});
set(ax2(2), 'Position', a2_pos{2});
set(ax2(3), 'Position', a2_pos{3});

% --------------------------------------------------------------------------
% Contact Area Schematic
%--------------------------------------------------------------------------

imshow('ContactAreaSchematic.png','Parent',ax3)
set(ax3, 'Position',[0.51 0.55 0.45 0.4]);

% --------------------------------------------------------------------------
% NP Capacity
%--------------------------------------------------------------------------

cla(ax4);

a4_pos = get(ax4, 'Position');

tcrs_per_cluster = [1,3,5,10,20];

col = [[0.1,0.1,0.5];[0.3,0.3,0.9];[0.3,0.5,1]];

hold(ax4,'on')
plotCapacity(ax4, tcrs_per_cluster, 8, 'np', col(1,:), true)
plotCapacity(ax4, tcrs_per_cluster, 14, 'np', col(2,:), true)
plotCapacity(ax4, tcrs_per_cluster, 20, 'np', col(3,:), true)
ylim(ax4,[0 310])
yticks(ax4, [0, 100, 200, 300])
ylabel(ax4, 'Capacity (number of NPs)');
xlabel(ax4, 'TCRs per Cluster (TPC)');

set(findall(ax4,'-property','FontSize'),'FontSize',10)

grid(ax4,'on');

set(ax4, 'Position', a4_pos,'box','off','tickdir','out');

l = legend(ax4,{'8 nm','14 nm','20 nm'},'box','off','Interpreter','Latex','Units','normalized','Position',[0.15 0.37 0.1 0.05],'FontSize',10,'FontWeight','bold');
title(l, "\textbf{NP Radius}",'Interpreter','latex','FontSize',10);

set(findall(ax4,'-property','Interpreter'),'Interpreter','Latex')

% --------------------------------------------------------------------------
% TCR Coverage
%--------------------------------------------------------------------------

cla(ax5(:));

a4_pos = get(ax5, 'Position');

np_radius = [8,14,20];
tcrs_per_cluster = [1,3,20];

col = [[0, 0, 0]; [0,0,1]; [1, 0, 0]];

for i = 1:length(np_radius)
    rNP = np_radius(i);
    title(ax5(i),['\textbf{',num2str(rNP),' nm}'],'Interpreter','latex','FontSize',14);

    for j = 1:length(tcrs_per_cluster)
        tpc = tcrs_per_cluster(j);
        load(['CapacityResults/Radius',num2str(rNP),'_TPC',num2str(tpc),'.mat'])
        hold(ax5(i),"on");
        plot(ax5(i),1:300, mean(Covered_TCRs,1), 'Color',col(j,:),'LineWidth',1.5);
    end
end
set(ax5(:),'box','off','tickdir','out')
ylim(ax5(:),[0 2.5]);
xlim(ax5(:), [0 300]);
xticks(ax5(:),[0 100 200 300]);
xticklabels(ax5(1:2),'');
xlabel(ax5(3), 'Bound NPs','FontSize',10);
ylabel(ax5(2), 'Covered TCRs per NP','FontSize',10);

%set(findall(ax5,'-property','FontSize'),'FontSize',10)
set(findall(ax5,'-property','Interpreter'),'Interpreter','Latex')

%--------------------------------------------------------------------------
% Time Series
%--------------------------------------------------------------------------

cla(ax6);

a6_pos = get(ax6, 'Position');

tcrs_per_cluster = [1,3,5,10,20];

col = [[0, 0, 0]; [0,0,1]; [1, 0, 0]];

hold(ax6,'on')

load("LongSims/r20/1TPC/v5/koff1/Uniform_rho100000.mat");
TCR_time_series = homo_bound_tcr{2};
plot(ax6, TCR_time_series(2,:), TCR_time_series(1,:),'Color',col(1,:));

load("LongSims/r20/3TPC/v5/koff1/Cluster_rho100000.mat");
TCR_time_series = cluster_bound_tcr{1};
plot(ax6, TCR_time_series(2,:), TCR_time_series(1,:),'Color',col(2,:));

clear("cluster_bound_tcr");
load("LongSims/r20/20TPC/v5/koff1/Cluster_rho100000.mat");
TCR_time_series = cluster_bound_tcr{2};
plot(ax6, TCR_time_series(2,:), TCR_time_series(1,:),'Color',col(3,:));

xlim(ax6, [1e2 1e4]);
xticklabels(ax6, {'2','4','6','8','10'})

xlabel(ax6,' Simulated Time ($10^3$ s)')
ylabel(ax6, 'Bound TCRs');

set(ax6,'tickdir','out','box','off');

set(findall(ax6,'-property','FontSize'),'FontSize',10)
set(findall(ax6,'-property','Interpreter'),'Interpreter','Latex')

% --------------------------------------------------------------------------
% Distributions Bound TCRs
%--------------------------------------------------------------------------

cla(ax7);

a7_pos = get(ax7, 'Position');

tcrs_per_cluster = [1,3,5,10,20];

col = [[0, 0, 0]; [0,0,1]; [1, 0, 0]];

hold(ax7,'on')

[h1,a1,v1] = TCR_Distributions(ax7, koff_vals(3), rho_vals(11), 5, 1, 'tcrs',col(1,:));
[h1,a1,v1] = TCR_Distributions(ax7, koff_vals(3), rho_vals(11), 5, 3, 'tcrs',col(2,:));
[h1,a1,v1] = TCR_Distributions(ax7, koff_vals(3), rho_vals(11), 5, 20, 'tcrs',col(3,:));

yticks(ax7,[0.02, 0.04, 0.06]);
ytickangle(ax7,0);
xlabel(ax7, 'Steady State Bound TCRs');
ylabel(ax7, 'Probability Density');

set(findall(ax7,'-property','FontSize'),'FontSize',10)
set(findall(ax7,'-property','Interpreter'),'Interpreter','Latex')

% Annotations
delete(findall(gcf,'type','annotation'));

annotation('textbox',[0.02 0.9 0.1 0.1],'string','(a)','FontSize',16, 'FontWeight','bold','Linestyle','none');
annotation('textbox',[0.48 0.9 0.1 0.1],'string','(b)','FontSize',16, 'FontWeight','bold','Linestyle','none');
annotation('textbox',[0.48 0.75 0.1 0.1],'string','(c)','FontSize',16, 'FontWeight','bold','Linestyle','none');

annotation('textbox',[0.02 0.52 0.1 0.1],'string','(d)','FontSize',16, 'FontWeight','bold','Linestyle','none');
annotation('textbox',[0.48 0.52 0.1 0.1],'string','(e)','FontSize',16, 'FontWeight','bold','Linestyle','none');

annotation('textbox',[0.02 0.22 0.1 0.1],'string','(f)','FontSize',16, 'FontWeight','bold','Linestyle','none');
annotation('textbox',[0.48 0.22 0.1 0.1],'string','(g)','FontSize',16, 'FontWeight','bold','Linestyle','none');

%% Save Figure
%print(f, '../Revised Figures/Figure1','-dpdf','-r0');