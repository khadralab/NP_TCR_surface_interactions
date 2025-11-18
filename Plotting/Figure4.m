%% Manuscript Figure 4
% Comparison of different NP valences
addpath('Plotting/Functions')
addpath('Functions/')
%%
% New Simulations 
koff_vals = [0.0001, 0.001, 0.01, 0.1];
kd_vals = koff_vals ./ 0.1;
rho_vals = linspace(-4,3.5,16);
bound_type = {'phos', 'bound', 'np', 'covered','fracphos', 'normBound', 'normPhos','rescaled'};

valence = [1,2,3,4,5,10];
tpc = [1, 3, 20];
np_radius = 20;


f = figure('Color','white','units','centimeters','Position',[55 1 20 24],...
    'PaperUnits','centimeters','PaperPosition', [1 1 20 24], 'Renderer','painters');

dr_left = [0.52 0.67 0.82];
dr_bot = [0.84, 0.70];
ec_left = [0.25, 0.5, 0.73];

for i=1:3
    ax1(i) = axes(f, 'Units','normalized','Position',[dr_left(i) 0.87 0.13 0.10]);
    ax2(i) = axes(f, 'Units','normalized','Position',[dr_left(i) 0.74 0.13 0.10]);
    ax3(i) = axes(f, 'Units','normalized','Position',[dr_left(i) 0.61 0.13 0.10]);
    ax5(i) = axes(f, 'Units','normalized','Position',[ec_left(i) 0.19 0.2 0.10]);
    ax6(i) = axes(f, 'Units','normalized','Position',[ec_left(i) 0.05 0.2 0.10]);
    if i==3
        ax4(1) = axes(f, 'Units','normalized','Position',[0.1 0.385 0.35 0.16]);
        ax4(2) = axes(f, 'Units','normalized','Position',[0.59 0.385 0.35 0.16]);
        ax7(1) = axes(f, 'Units','normalized','Position',[0.05 0.20 0.13 0.07]);
        ax7(2) = axes(f, 'Units','normalized','Position',[0.05 0.1 0.13 0.07]);
    end
end
% Color and Linestyle
col = gray(length(kd_vals)+2);
hol = winter(length(kd_vals)+2);
sol = autumn(length(kd_vals)+2);

linestyle = {'-','--'};
markerstyle = {'s','d'};
colorstyle = {col, hol, sol};

%--------------------------------------------------------------------------
% Dose Response Curves
%--------------------------------------------------------------------------

a1_pos = get(ax1, 'Position');
a2_pos = get(ax2, 'Position');
a3_pos = get(ax3, 'Position');

kp1 = [0.28];
%kp1 = [0.0001,2;0.001,2; 0.01,2;0.1,2];
%kp2 = [0.01,1; 0.01,2; 0.01,3; 0.01,4];
surfaces = [1,3,20];

for i=1:length(kd_vals)
    for j=1:length(surfaces)
        %DoseResponsePlots(ax2(i), kd_vals(i), kp1(j,:), rho_vals, '33x3',...
        %5, bound_type{7}, colorstyle{2}(j,:), false, false, np_radius);
        %hold(ax2(i),"on");

        %DoseResponsePlots(ax3(i), kd_vals(i), kp2(j,:), rho_vals, '10x20',...
        %5, bound_type{7}, colorstyle{2}(j,:), false, false, np_radius);
        %hold(ax3(i),"on");
        DoseResponsePlots(ax1(j), kd_vals(i), [kp1,5], rho_vals, surfaces(j), 5, bound_type{6}, colorstyle{j}(i,:), false, false, np_radius);
        hold(ax1(j),"on");

        DR_KPR_negativeFeedback(ax2(j), kd_vals(i), kp1, rho_vals, surfaces(j), 5, colorstyle{j}(i,:), false, false, np_radius, false);
        hold(ax2(j),"on");
        DR_KPR_negativeFeedback(ax3(j), kd_vals(i), kp1, rho_vals, surfaces(j), 5, colorstyle{j}(i,:), false, false, np_radius, true);
        hold(ax3(j),"on");
    end
end

for i=1:3
    set(ax1(i), 'Position', a1_pos{i});
    set(ax2(i), 'Position', a2_pos{i});
    set(ax3(i), 'Position', a3_pos{i});
    title(ax1(i),['TPC = ',num2str(surfaces(i))]);
end

set([ax1(:) ax2(:) ax3(:)],'xscale','log','tickdir','out','box','off')
ylim([ax1(:) ax2(:) ax3(:)], [0 300]);
xlim([ax1(:) ax2(:) ax3(:)], [0.0001 15000]);
xlabel([ax3(2)], 'pMHC Concentration (log$_{10}$ molecules/second)');
ylabel([ax2(1) ax3(1)], 'Phos TCRs')
ylabel([ax1(1)], 'Bound TCRs')


set([ax1(:) ax2(:) ax3(:)],'GridLineWidth',0.5,'Linewidth',2)
xlbl = 10.^(-3:2:3);
xlabels = {'-3', '-1', '1', '3'};
xticks([ax1(:) ax2(:) ax3(:)],xlbl);
xticklabels([ax1(:) ax2(:) ax3(:)], xlabels)
%xsecondarylabel([ax(2) ax(5)], '')
set([ax1(:) ax2(:) ax3(:)],'TickLabelInterpreter', 'Latex')
yticklabels([ax1(2:3) ax2(2:3) ax3(2:3)], '')

set(findall([ax1 ax2 ax3],'-property','FontSize'),'FontSize',10)
set(findall([ax1 ax2 ax3],'-property','Interpreter'),'Interpreter','Latex')

% -------------------------------------------------------------------------
% EMax versus EC50 (Row 3)
%--------------------------------------------------------------------------

cla(ax4(1));
a4_pos = get(ax4, 'Position');

koff_vals = [0.0001, 0.001, 0.01, 0.1];
kd_vals = koff_vals ./ 0.1;
kp1 = [0.001,3; 0.005,3; 0.01,3; 1,0;];

surfaces = [1,3, 5, 10, 20];
surfaces2 = {1, '5x20', '10x20',20};

colors = parula(length(koff_vals)+2);

line1 = Ec50_Emax_KPR(ax4(1), kd_vals, rho_vals, 3, 5, kp1, 'phos', col);

set(ax4(1), 'Position', a4_pos{1}+[0, -0.02, 0, 0]);

xlabel([ax4(1)],'EC50');
ylabel([ax4(1)],'EMax');
xlim([ax4(1)], [10^-2.5 10^3]);
ylim([ax4(1)], [0 300]);
xticks(ax4(1), [10^-2 10^0 10^2]);
set([ax4(1)], 'xscale','log','yscale','linear');
grid([ax4(1)],'on');

l=legend(ax4(1),{'$10^{-4}$','$10^{-3}$','$10^{-2}$','$10^{-1}$'},'Interpreter','latex','Box','off',...
    'Orientation','vertical','Direction','reverse','NumColumns',2,'Position',[0.38 0.48 0.03 0.03]);
title(l,"$k_{off}$",'Interpreter','latex');

set(findall(ax4(1),'-property','FontSize'),'FontSize',10)
set(findall(ax4(1),'-property','Interpreter'),'Interpreter','Latex')

%--------------------------------------------------------------------------
% Discriminatory Power (Row 3)
%--------------------------------------------------------------------------

cla(ax4(2));
a4_pos = get(ax4, 'Position');

kpn = [0;1;3;5];
kp2 = 0.01*ones(4,1);
kp2 = [kp2 ,kpn];

plot_DP_KPR(ax4(2), 1, 5, 3, kp2);

set(ax4(2), 'Position', a4_pos{2}+[0,-0.02,0,0]);

xlabel([ax4(2)],'$K_D$');
ylabel([ax4(2)],'pMHC Cncentration');
xlim([ax4(2)], [10^-3.1 10^0.1]);
ylim([ax4(2)], [10^-4 10^2.2]);
yticks(ax4(2), [10^-4 10^-2 10^0 10^2]);
xticks(ax4(2), [10^-3 10^-2 10^-1 10^0]);
set([ax4(2)], 'xscale','log','yscale','log');
grid([ax4(2)],'on');

l=legend(ax4(2),{'$0$','$1$','$3$','$5$'},'Interpreter','latex','Box','off',...
    'Orientation','vertical','Direction','reverse','NumColumns',2,'Position',[0.83 0.38 0.03 0.03]);
title(l,"$N$",'Interpreter','latex');

set(findall(ax4(2),'-property','FontSize'),'FontSize',10)
set(findall(ax4(2),'-property','Interpreter'),'Interpreter','Latex')

% -----------------------------------------------------------------------
% Activation Heatmaps
%-------------------------------------------------------------------------

valence_vals = [5];
surfaces = [1,3,20];
threshold_array = [15];

cla(ax5(:)); cla(ax6(:));

a6_pos = get(ax6, 'Position');
a5_pos = get(ax5, 'Position');

col = hot(20);
col = flipud(col);

% Uniform (1 TPC) Surface
KPR_heatmap(ax5(1), 10, 1, 0.01, false, false);
KPR_heatmap(ax5(2), 10, 1, 0.01, true, false);
KPR_heatmap(ax5(3), 10, 1, 0.01, true, true);

% Clustered (20 TPC) Surface
KPR_heatmap(ax6(1), 10, 20, 0.01, false, false);
KPR_heatmap(ax6(2), 10, 20, 0.01, true, false);
KPR_heatmap(ax6(3), 10, 20, 0.01, true, true);

set([ax5(:) ax6(:)],'xscale','log','yscale','log','tickdir','out','box','off')
grid([ax5(:) ax6(:)], 'on');

ylabel([ax5(1) ax6(1)], '$K_D$');
xlabel(ax6(2), 'pMHC Concentration (molecules/second)');
yticks([ax5(:) ax6(:)], [10^-4 10^-3 10^-2 10^-1]);
xticks([ax5(:) ax6(:)], [10^-3 10^-1 10^1 10^3]);
xticklabels(ax5(:), '');
yticklabels([ax5(2:3) ax6(2:3)],'');
title(ax5(1), 'No Proofreading');
title(ax5(2), 'Canonical KPR');
title(ax5(3), 'Modified KPR');

%cb1 = colorbar(ax5(3), 'Ticks',[0, 1, 5, 10, 50, 100],'Position',[0.9 0.15 0.01 0.05],'HandleVisibility','on');
%cb1.Label.String = 'Activation';

for i=1:3
    set(ax5(i), 'Position', a5_pos{i});
    set(ax6(i), 'Position', a6_pos{i});
    colormap(ax5(i), col);
    colormap(ax6(i), col);
end

set([ax5(:) ax6(:)],'clim',[0.1 100]);
set([ax5(:) ax6(:)],'ColorScale','log');

set(findall([ax5 ax6],'-property','FontSize'),'FontSize',10)
set(findall([ax5 ax6],'-property','Interpreter'),'Interpreter','Latex')


%% --------------------------------------------------------------------------
% Surfaces (1 TPC & 20 TPC)
%--------------------------------------------------------------------------
a7_pos = get(ax7, 'Position');

rSurf = 1000;
tcr_per_cluster = 1;
num_tcr = 300;
max_clusters = floor(num_tcr/tcr_per_cluster);
num_clusters = [5,10,20];
cluster_radius = 50;
rTCR = 5;
rNP = 20;

col = [[0, 0, 0]; [1, 0, 0]];

tcr_pos = gen_clusters(rSurf, cluster_radius, 20, 0, num_tcr);
PlotSurface(ax7(1), rSurf, tcr_pos, rTCR, [], rNP, col(1,:));

tcr_pos = gen_clusters(rSurf, cluster_radius, 15, 20, num_tcr);
PlotSurface(ax7(2), rSurf, tcr_pos, rTCR, [], rNP, col(2,:));

set(ax7(1), 'Position', a7_pos{1}+[-0.02 0 0 0.02]);
set(ax7(2), 'Position', a7_pos{2}+[-0.02 -0.04 0 0.02]);

set(findall(ax7(:),'-property','FontSize'),'FontSize',10)
set(findall(ax7(:),'-property','Interpreter'),'Interpreter','Latex')

%% Annotations (Panel labels)
delete(findall(gcf,'type','annotation'));

annotation('textbox',[0.02 0.9 0.1 0.1],'string','(a)','FontSize',16, 'FontWeight','bold','Linestyle','none');
annotation('textbox',[0.02 0.7 0.1 0.1],'string','(b)','FontSize',16, 'FontWeight','bold','Linestyle','none');
annotation('textbox',[0.4 0.9 0.1 0.1],'string','(c)','FontSize',16, 'FontWeight','bold','Linestyle','none');

annotation('textbox',[0.02 0.47 0.1 0.1],'string','(d)','FontSize',16, 'FontWeight','bold','Linestyle','none');
annotation('textbox',[0.47 0.47 0.1 0.1],'string','(e)','FontSize',16, 'FontWeight','bold','Linestyle','none');

annotation('textbox',[0.02 0.23 0.1 0.1],'string','(f)','FontSize',16, 'FontWeight','bold','Linestyle','none');

%% Save Figure
print(f, '../Revised Figures/Figure4_noSchem','-dpdf','-r0');