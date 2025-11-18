%% Manuscript Figure 2
% 1) Comparison of different surfaces
addpath('Plotting/Functions')
addpath('Functions/')

%%
f = figure('Color','white','units','centimeters','Position',[55 1 20 24],...
    'PaperUnits','centimeters','PaperPosition', [1 1 20 24], 'Renderer','painters');

rhs_left = [0.52 0.67 0.82];
lhs_left = [0.05, 0.2, 0.35];
mi_left = [0.52 0.65 0.78];
dr_bot = [0.84, 0.70];
ec_left = [0.15, 0.4, 0.65];

for i=1:3
    ax2(i) = axes(f, 'Units','normalized','Position',[rhs_left(i) 0.84 0.13 0.10]); % DR Curves
    ax7(i) = axes(f, 'Units','normalized','Position',[mi_left(i) 0.6 0.12 0.15]); % MI Heatmaps
    ax8(i) = axes(f, 'Units','normalized','Position',[lhs_left(i)+0.02 0.18 0.13 0.10]); % Surfaces Mixtures
    ax9(i) = axes(f, 'Units','normalized','Position',[lhs_left(i)+0.02 0.05 0.13 0.10]); % DR Curves Mixtures

    if i==3
        ax3 = axes(f, 'Units','normalized','Position',[0.1 0.37 0.3 0.15]); % Threshold conc. versus Kd
        ax4 = axes(f, 'Units','normalized','Position',[0.52 0.37 0.20 0.15]); % Low affinity DP
        ax5 = axes(f, 'Units','normalized','Position',[0.75 0.37 0.20 0.15]); % High affinity DP
        ax6 = axes(f, 'Units','normalized','Position',[0.1 0.6 0.3 0.15]); % EMax v. EC50
        ax10 = axes(f, 'Units','normalized','Position',[0.6 0.05 0.3 0.22]); % EMax v. EC50
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


% --------------------------------------------------------------------------
% Dose Response Curves
%--------------------------------------------------------------------------
a2_pos = get(ax2, 'Position');

koff_vals = [0.0001, 0.001, 0.01];
kd_vals = koff_vals ./ 0.1;

for i=1:length(kd_vals)
    for j=1:length(valence)
        for k=1:length(tpc)
            utop{j,k} = DoseResponsePlots(ax2(i), kd_vals(i), kp, rho_vals, tpc(k),...
            valence(j), bound_type{6}, colorstyle{k}(j,:), false, false, np_radius);
            hold(ax2(i),"on");
        end
    end
end

for i=1:3
    set(ax2(i), 'Position', a2_pos{i});
end

% Format Axes

set(ax2(:),'xscale','log','tickdir','out','box','off')
ylim(ax2(:), [0 300]);
xlim(ax2(:), [0.0001 15000]);
xlabel(ax2(2), 'pMHC Concentration (log$_{10}$ molecules/second)');
ylabel(ax2(1), 'Bound TCRs')

set(ax2(:),'GridLineWidth',0.5,'Linewidth',2)
xlbl = 10.^(-3:2:3);
xlabels = {'-3', '-1', '1', '3'};
xticks(ax2(:),xlbl);
xticklabels(ax2(:), xlabels)
%xsecondarylabel([ax(2) ax(5)], '')
set(ax2(:),'TickLabelInterpreter', 'Latex')
yticklabels(ax2(2:3), '')

for i=1:3
     title(ax2(i),['$k_{off} = 10^{',num2str(log10(koff_vals(i))),'}$']);
end

set(findall(ax2,'-property','FontSize'),'FontSize',10)
set(findall(ax2,'-property','Interpreter'),'Interpreter','Latex')
%--------------------------------------------------------------------------
% Discriminatory Power
%--------------------------------------------------------------------------
cla([ax3, ax4, ax5, ax6]);

valence_vals = [5];
surfaces = [1,3,5,10,20];
threshold_array = [10];

a3_pos = get(ax3, 'Position');
a4_pos = get(ax4, 'Position');
a5_pos = get(ax5, 'Position');

plot_Discriminatory_Power(ax3, ax4, ax5, threshold_array, valence_vals, surfaces, kp);

set(ax3, 'Position', a3_pos);
set(ax4, 'Position', a4_pos);
set(ax5, 'Position', a5_pos);

%ylim(ax3, [10^-4 10^3]);
%xlim(ax3, [10^-3.7 10^0.5]);

ylim([ax4 ax5], [0 2]);
xlim([ax4 ax5], [0 20]);
xticks([ax4 ax5],[1 3 5 10 20]);
yticks([ax4 ax5], linspace(0,3,4));

xlabel(ax3, '$K_D$');
ylabel(ax3, 'pMHC Concentration');

xlabel([ax4 ax5], 'TPC');
ylabel(ax4, 'Discriminatory Power');
yticklabels(ax5, '');
box(ax3, 'off');

set(findall([ax3 ax4 ax5],'-property','FontSize'),'FontSize',10)
set(findall([ax3, ax4, ax5],'-property','Interpreter'),'Interpreter','Latex')
%--------------------------------------------------------------------------
% EMax versus EC50
%--------------------------------------------------------------------------
%koff_vals = [0.0001, 0.0003, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3];
koff_vals = [0.0001, 0.001, 0.01, 0.1];
kd_vals = koff_vals ./ 0.1;

surfaces = [1,3,5,10,20];
colors = parula(length(koff_vals)+2);

a6_pos = get(ax6, 'Position');

line1 = Ec50_Emax(ax6, kd_vals, rho_vals, surfaces, valence, kp, 'bound', colors);
%line2 = Ec50_Emax(ax6, kd_vals, rho_vals, 3, valence, 'bound', hol);
%line3 = Ec50_Emax(ax6, kd_vals, rho_vals, 20, valence, 'bound', sol);

set(ax6,'Position', a6_pos);

xlabel(ax6,'EC50');
ylabel(ax6,'EMax');
xlim(ax6, [10^-2.5 10^2]);
ylim(ax6, [0 300]);
set(ax6, 'xscale','log','yscale','linear');
grid(ax6,'on');

l=legend(ax6,{'$10^{-4}$','$10^{-3}$','$10^{-2}$','$10^{-1}$'},'Interpreter','latex','Box','off',...
    'Orientation','vertical','Direction','reverse','NumColumns',2,'Position',[0.17 0.625 0.03 0.03]);
title(l,"$k_{off}$",'Interpreter','latex');

set(findall(ax6,'-property','FontSize'),'FontSize',10)
set(findall(ax6,'-property','Interpreter'),'Interpreter','Latex')

%--------------------------------------------------------------------------
% MI heatmap
%--------------------------------------------------------------------------
a7_pos = get(ax7,'Position');
surfaces = [1,3,20];
i=1;
normtrue = false;
for t = 1:length(surfaces)
    TPC = surfaces(t);

    plot_MI_heatmap(ax7(t), TPC, "kd", false, "mi");

    set(ax7(t),'Position',a7_pos{t})
    
    %plot_MI_heatmap(ax(t+4), TPC, "kd", false, "mi");
end

yticklabels(ax7(2:3), '');
xtickangle(ax7(:), 0);
xlabel(ax7(2), 'NP Concentration (log$_{10}$ molecules/second)')
ylabel(ax7(1) ,'NP valence');

set(findall(ax7,'-property','FontSize'),'FontSize',10)
set(findall(ax7,'-property','Interpreter'),'Interpreter','Latex')

%--------------------------------------------------------------------------
% Plot Surfaces
%--------------------------------------------------------------------------
rSurf = 1000;
tcr_per_cluster = 20;
num_tcr = 300;
max_clusters = floor(num_tcr/tcr_per_cluster);
num_clusters = [5,10,20];
cluster_radius = 50;
rTCR = 5;
rNP = 20;

col = [[0, 0, 0]; [1, 0.4, 0]; [1, 0, 0.4]; [1, 0, 0]];

tcr_pos = gen_clusters(rSurf, cluster_radius, 20, 0, num_tcr);
PlotSurface(ax8(1), rSurf, tcr_pos, rTCR, [], rNP, col(1,:));

for i=2:3
    tcr_pos = gen_clusters(rSurf, cluster_radius, num_clusters(i-1), tcr_per_cluster, num_tcr);
    PlotSurface(ax8(i), rSurf, tcr_pos, rTCR, [], rNP, col(i,:));
end

set(findall(ax8(:),'-property','FontSize'),'FontSize',16)
set(findall(ax8(:),'-property','Interpreter'),'Interpreter','Latex')

% --------------------------------------------------------------------------
% DR Mixtures
%--------------------------------------------------------------------------

a9_pos = get(ax9, 'Position');

koff_vals = [0.0001, 0.001, 0.01];
kd_vals = koff_vals ./ 0.1;

for i=1:length(kd_vals)
    DoseResponsePlots(ax9(i), kd_vals(i), kp, rho_vals, 1,...
    5, bound_type{7}, col(1,:), false, false, np_radius);
    hold(ax9(i),"on");

    DoseResponsePlots(ax9(i), kd_vals(i), kp, rho_vals, '5x20',...
    5, bound_type{7}, col(2,:), false, false, np_radius);

    DoseResponsePlots(ax9(i), kd_vals(i), kp, rho_vals, '10x20',...
    5, bound_type{7}, col(3,:), false, false, np_radius);

    DoseResponsePlots(ax9(i), kd_vals(i), kp, rho_vals, 20,...
    5, bound_type{7}, col(4,:), false, false, np_radius);
end

for i=1:3
    set(ax9(i), 'Position', a9_pos{i});
end

% Format Axes

set(ax9(:),'xscale','log','tickdir','out','box','off')
ylim(ax9(:), [0 300]);
xlim(ax9(:), [0.0001 15000]);
xlabel(ax9(2), 'pMHC Concentration (log$_{10}$ molecules/second)');
ylabel(ax9(1), 'Bound TCRs')

set(ax9(:),'GridLineWidth',0.5,'Linewidth',2)
xlbl = 10.^(-3:2:3);
xlabels = {'-3', '-1', '1', '3'};
xticks(ax9(:),xlbl);
xticklabels(ax9(:), xlabels)
%xsecondarylabel([ax(2) ax(5)], '')
set(ax9(:),'TickLabelInterpreter', 'Latex')
yticklabels(ax9(2:3), '')

for i=1:3
     title(ax9(i),['$k_{off} = 10^{',num2str(log10(koff_vals(i))),'}$']);
end

set(findall(ax9,'-property','FontSize'),'FontSize',10)
set(findall(ax9,'-property','Interpreter'),'Interpreter','Latex')

% --------------------------------------------------------------------------
% EC50 - EMax Mixtures
%--------------------------------------------------------------------------

a10_pos = get(ax10, 'Position');

koff_vals = [0.0001, 0.001, 0.01, 0.1];
kd_vals = koff_vals ./ 0.1;

surfaces2 = {1, '5x20', '10x20',20};
colors = parula(length(koff_vals)+2);

line1 = Ec50_Emax_Mixtures(ax10, kd_vals, rho_vals, surfaces2, 5, kp, 'bound', col);
hold(ax10,'on');

l=legend(ax10,{'$10^{-4}$','$10^{-3}$','$10^{-2}$','$10^{-1}$'},'Interpreter','latex','Box','off',...
    'Orientation','vertical','Direction','reverse','NumColumns',2,'Position',[0.83 0.23 0.03 0.03]);
title(l,"$k_{off}$",'Interpreter','latex');

set(ax10, 'Position', a10_pos);

xlabel(ax10,'EC50');
ylabel(ax10,'EMax');
xlim(ax10, [10^-2.5 10^2]);
ylim(ax10, [0 300]);
set(ax10, 'xscale','log','yscale','linear');
grid(ax10,'on');

set(findall(ax10,'-property','FontSize'),'FontSize',10)
set(findall(ax10,'-property','Interpreter'),'Interpreter','Latex')

% --------------------------------------------------------------------------
% Annotations (panel labels)
%--------------------------------------------------------------------------
delete(findall(gcf,'type','annotation'));

annotation('textbox',[0.01 0.9 0.1 0.1],'string','(a)','FontSize',14, 'FontWeight','bold','Linestyle','none');
annotation('textbox',[0.43 0.9 0.1 0.1],'string','(b)','FontSize',14, 'FontWeight','bold','Linestyle','none');
annotation('textbox',[0.01 0.69 0.1 0.1],'string','(c)','FontSize',14, 'FontWeight','bold','Linestyle','none');
annotation('textbox',[0.45 0.69 0.1 0.1],'string','(d)','FontSize',14, 'FontWeight','bold','Linestyle','none');
annotation('textbox',[0.01 0.46 0.1 0.1],'string','(e)','FontSize',14, 'FontWeight','bold','Linestyle','none');
annotation('textbox',[0.45 0.46 0.1 0.1],'string','(f)','FontSize',14, 'FontWeight','bold','Linestyle','none');
annotation('textbox',[0.72 0.46 0.1 0.1],'string','(g)','FontSize',14, 'FontWeight','bold','Linestyle','none');
annotation('textbox',[0.01 0.22 0.1 0.1],'string','(h)','FontSize',14, 'FontWeight','bold','Linestyle','none');
annotation('textbox',[0.51 0.22 0.1 0.1],'string','(i)','FontSize',14, 'FontWeight','bold','Linestyle','none');

%% Save Figure
%print(f, '../Revised Figures/Figure2-test','-dpdf','-r0');
