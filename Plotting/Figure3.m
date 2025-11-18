%% Manuscript Figure 3
% Comparison of different NP valences
addpath('Plotting/Functions')
addpath('Functions/')
%%

f1 = figure('Color','white','units','centimeters','Position',[20 1 20 24],...
    'PaperUnits','centimeters','PaperPosition', [1 1 20 24], 'Renderer','painters');

dr_left = [0.17 0.32 0.53];
dr_bot = [0.84, 0.70];
ec_left = [0.17, 0.45, 0.72];

for i=1:3
    ax2(i) = axes(f1, 'Units','normalized','Position',[dr_left(i) dr_bot(1) 0.13 0.10]); % DR 1TPC
    ax3(i) = axes(f1, 'Units','normalized','Position',[dr_left(i) dr_bot(2) 0.13 0.10]); % DR 20 TPC
    ax4(i) = axes(f1, 'Units','normalized','Position',[ec_left(i) 0.47 0.22 0.12]); % EC50-EMax
    
    if i==3
        ax1(1) = axes(f1, 'Units','normalized','Position',[0.02 dr_bot(1)+0.01 0.10 0.08]); % Surfaces
        ax1(2) = axes(f1, 'Units','normalized','Position',[0.02 dr_bot(2)+0.01 0.10 0.08]);
        ax5(1) = axes(f1, 'Units','normalized','Position',[0.73 dr_bot(1) 0.15 0.10]); % MI 
        ax5(2) = axes(f1, 'Units','normalized','Position',[0.73 dr_bot(2) 0.15 0.10]); % MI 
        ax6 = axes(f1, 'Units','normalized','Position',[0.17 0.12 0.45 0.27]); % DP
        ax7 = axes(f1, 'Units','normalized','Position',[0.72 0.12 0.22 0.1]);
        ax8 = axes(f1, 'Units','normalized','Position',[0.72 0.29 0.22 0.1]);
        break
    end
end

% New Simulations
koff_vals = [0.0001, 0.001, 0.01];
kd_vals = koff_vals ./ 0.1;
rho_vals = linspace(-4,3.5,16);
bound_type = {'phos', 'bound', 'np', 'covered','fracphos', 'normBound', 'normPhos','rescaled'};

valence = [1,2,3,4,5,10];
tpc = [1, 3, 20];
kp = [0.01,2];
np_radius = 20;

% Color and Linestyle
col = gray(length(valence)+2);
hol = winter(length(valence)+2);
sol = autumn(length(valence)+2);

linestyle = {'-','--'};
markerstyle = {'s','d'};
colorstyle = {col, hol, sol};
%%
%--------------------------------------------------------------------------
% Surfaces (1 TPC & 20 TPC)
%--------------------------------------------------------------------------
a1_pos = get(ax1, 'Position');

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
PlotSurface(ax1(1), rSurf, tcr_pos, rTCR, [], rNP, col(1,:));

tcr_pos = gen_clusters(rSurf, cluster_radius, 15, 20, num_tcr);
PlotSurface(ax1(2), rSurf, tcr_pos, rTCR, [], rNP, col(2,:));

set(findall(ax1(:),'-property','FontSize'),'FontSize',16)
set(findall(ax1(:),'-property','Interpreter'),'Interpreter','Latex')
%--------------------------------------------------------------------------
% Dose Response Curves
%--------------------------------------------------------------------------

a2_pos = get(ax2, 'Position');
a3_pos = get(ax3, 'Position');

for i=1:2
    hold(ax2(i),"on");
    hold(ax3(i),"on");
    for j=1:length(valence)
        DoseResponsePlots(ax2(i), kd_vals(i), kp, rho_vals, tpc(1),...
        valence(j), bound_type{6}, colorstyle{1}(j,:), false, false, np_radius);
        
        DoseResponsePlots(ax3(i), kd_vals(i), kp, rho_vals, tpc(3),...
        valence(j), bound_type{6}, colorstyle{3}(j,:), false, false, np_radius);
    end
    title(ax2(i),['$K_{D} = ',num2str(kd_vals(i)),'$']);
end

hold([ax2(3) ax3(3)],"on")
for j=1:length(valence)
    DoseResponsePlots(ax2(3), kd_vals(1), kp, rho_vals, tpc(1),...
    valence(j), bound_type{3}, colorstyle{1}(j,:), false, false, np_radius);
    
    DoseResponsePlots(ax3(3), kd_vals(1), kp, rho_vals, tpc(3),...
    valence(j), bound_type{3}, colorstyle{3}(j,:), false, false, np_radius);
end
title(ax2(3),['$K_{D} = ',num2str(kd_vals(1)),'$']);

% Format Axes
for i=1:3
    set(ax2(i), 'Position', a2_pos{i});
    set(ax3(i), 'Position', a3_pos{i});
end

set([ax2(:) ax3(:)],'xscale','log','tickdir','out','box','off')
ylim([ax2(:) ax3(:)], [0 300]);
xlim([ax2(:) ax3(:)], [0.0001 15000]);
xlabel([ax3(1)], 'pMHC Concentration (log$_{10}$)', 'Position',[2e5 -80]);
xlabel([ax3(3)], 'NP Concentration (log$_{10}$)');
ylabel([ax3(1)], 'Bound TCRs','Units','centimeters','Position', [-0.8 3]);
ylabel([ax3(3)], 'Bound NPs', 'Units','centimeters','Position', [-0.2 3]);

set([ax2(:) ax3(:)],'GridLineWidth',0.5,'Linewidth',2)
xlbl = 10.^(-3:2:3);
xlabels = {'-3', '-1', '1', '3'};
xticks([ax2(:) ax3(:)],xlbl);
xticklabels([ax2(:) ax3(:)], xlabels)
%xsecondarylabel([ax(2) ax(5)], '')
set([ax2(:) ax3(:)],'TickLabelInterpreter', 'Latex')
yticklabels([ax2(2:3) ax3(2:3)], '')

set(findall([ax2 ax3],'-property','FontSize'),'FontSize',10)
set(findall([ax2 ax3],'-property','Interpreter'),'Interpreter','Latex')

set(f1, 'Units', 'centimeters')
pos = get(f1, 'Position');
set(f1,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])

%--------------------------------------------------------------------------
% EMax versus EC50
%--------------------------------------------------------------------------
cla(ax4(:));
a4_pos = get(ax4, 'Position');

koff_vals = [0.0001, 0.001, 0.01, 0.1];
kd_vals = koff_vals ./ 0.1;

colors = parula(length(koff_vals)+2);

line1 = Ec50_Emax_valence(ax4(1), kd_vals, rho_vals, 1, valence, 'bound', col);
hold(ax4(1),"on");
line2 = Ec50_Emax_valence(ax4(2), kd_vals, rho_vals, 3, valence, 'bound', sol);
hold(ax4(2),"on");
line2 = Ec50_Emax_valence(ax4(3), kd_vals, rho_vals, 20, valence, 'bound', sol);
hold(ax4(3),"on");

set(findall(ax4(:),'-property','FontSize'),'FontSize',10)

l=legend(ax4(1),{'$10^{-4}$','$10^{-3}$','$10^{-2}$','$10^{-1}$'},'Interpreter','latex','box','off','NumColumns',2,...
    'Position',[0.335 0.555 0.01 0.01]);
title(l,"$k_{off}$",'Interpreter','latex','FontSize',12);

surfaces = [1,3,20];
for i=1:3
    set(ax4(i), 'Position', a4_pos{i});
    title(ax4(i), [num2str(surfaces(i)),' TPC'],'FontSize',12);
end

xlabel(ax4(2),'EC50 (NP Concentration)');
ylabel(ax4(1),'EMax (Bound TCRs)');
yticklabels(ax4(2:3), '');
xlim(ax4(:), [10^-2.5 10^2]);
xticks(ax4(:), [10^-2, 10^0, 10^2]);
ylim(ax4(:), [0 300]);
set(ax4(:), 'xscale','log','yscale','linear');
grid(ax4(:),'on');

set(findall(ax4(:),'-property','Interpreter'),'Interpreter','Latex')

%--------------------------------------------------------------------------
% Heatmap
%--------------------------------------------------------------------------
cla(ax5(:));
a5_pos = get(ax5, 'Position');

surfaces = [1,20];
i=1;
normtrue = false;
for t = 1:length(surfaces)
    TPC = surfaces(t);

    plot_MI_heatmap(ax5(t), TPC, "valence", false, "norm");
    
    %plot_MI_heatmap(ax(t+4), TPC, "kd", false, "mi");
end

for i=1:2
    set(ax5(i), 'Position', a5_pos{i});
end

set(findall(ax5(:),'-property','FontSize'),'FontSize',8)

cb = colorbar(ax5(2), 'Ticks',[0,0.05,0.1,0.15,0.2,0.25],'Units','centimeters','Position',[18.5 18 0.3 3]);
cb.Label.String = 'Mutual Information (bits)';
cb.Label.FontSize = 10;

xticklabels(ax5(1), '');
ylabel(ax5(1), '$K_D$','Units','centimeters','Position',[-0.7 -0.5],'FontSize',10);
xlabel(ax5(2), 'NP Concentration','FontSize',10);
title(ax5(1), '1 TPC','FontSize',12);
title(ax5(2), '20 TPC','FontSize',12);

set(findall(ax5(:),'-property','Interpreter'),'Interpreter','Latex')

%--------------------------------------------------------------------------
% Discriminatory Power
%--------------------------------------------------------------------------
a6_pos = get(ax6, 'Position');
a7_pos = get(ax7, 'Position');
a8_pos = get(ax8, 'Position');

valence_vals = [2,3,4,5,6,10];
threshold_array = [10];

plot_DP_valence(ax6, ax7, ax8, threshold_array, valence_vals, 20);

set(ax6, 'Position', a6_pos);
set(ax7, 'Position', a7_pos);
set(ax8, 'Position', a8_pos);

xlim(ax6, [10^-3.2 10^0.2]);
ylim(ax6, [10^-4 10^2.5]);
xlabel(ax6, '$K_D$','Interpreter','latex');
ylabel(ax6, 'pMHC Concentration');
l = legend(ax6, {'2', '5', '10'},'Location','northwest','box','off','Orientation','horizontal','FontSize',12);
title(l,'$v$','Interpreter','latex','FontSize',16);

ylim([ax7 ax8], [0 2.5]);
xlim([ax7 ax8], [1 10]);
xticks([ax7 ax8],1:10);
yticks([ax7 ax8], linspace(0,2,3));
xticklabels([ax7 ax8],{'','2','','4','','6','','','','10'})
xtickangle([ax7 ax8],0);
%xticklabels(ax8,'')
xlabel(ax7, 'NP Valence ($v$)');
ylabel(ax7, 'Discriminatory Power','units','centimeters','Position',[-0.7 3]);

set(findall([ax6 ax7 ax8],'-property','FontSize'),'FontSize',10)
set(findall([ax6, ax7, ax8],'-property','Interpreter'),'Interpreter','Latex')

%% Annotations (Panel labels)
delete(findall(gcf,'type','annotation'));

annotation('textbox',[0.04 0.88 0.1 0.1],'string','(a)','FontSize',16, 'FontWeight','bold','Linestyle','none');
annotation('textbox',[0.47 0.88 0.1 0.1],'string','(b)','FontSize',16, 'FontWeight','bold','Linestyle','none');
annotation('textbox',[0.67 0.88 0.1 0.1],'string','(c)','FontSize',16, 'FontWeight','bold','Linestyle','none');

annotation('textbox',[0.04 0.53 0.1 0.1],'string','(d)','FontSize',16, 'FontWeight','bold','Linestyle','none');
annotation('textbox',[0.4 0.53 0.1 0.1],'string','(e)','FontSize',16, 'FontWeight','bold','Linestyle','none');
annotation('textbox',[0.67 0.53 0.1 0.1],'string','(f)','FontSize',16, 'FontWeight','bold','Linestyle','none');

annotation('textbox',[0.04 0.33 0.1 0.1],'string','(g)','FontSize',16, 'FontWeight','bold','Linestyle','none');
annotation('textbox',[0.67 0.33 0.1 0.1],'string','(h)','FontSize',16, 'FontWeight','bold','Linestyle','none');

%% Save Figure
print(f1, '../Revised Figures/Figure3','-dpdf','-r0');