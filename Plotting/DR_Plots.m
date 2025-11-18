%% Dose-Response Curves
% Plots of results from DoseResponse.m for: (r,v, koff) = (20,5,0.05)
clear all
addpath('Plotting/Functions')

bound_type = {'phos', 'bound', 'np', 'covered','fracphos', 'xPMHC', 'rescaled'};
btype = 7;                          % 1 for phos; 2 for bound; 3 for NPs
tcrs_per_cluster = 1;
tcrs_per_cluster2 = 3;
tcrs_per_cluster3 = 5;
tcrs_per_cluster4 = 20;
valence = 5;

% Original Simulations
koff_vals = [0.01, 0.02, 0.03, 0.05];
kd_vals = koff_vals ./ 0.1;
rho_vals = linspace(-2,5,15);           % Full range of NP-concentrations.
rho_vals = rho_vals(1:12);              % Select the x-axis concentrations to illustrate


% New Simulations
koff_vals = [0.0001, 0.001, 0.01, 0.1];
kd_vals = koff_vals ./ 0.1;
rho_vals = linspace(-4,3.5,16);

%Plotting
hols = autumn(length(koff_vals)+1);
hols = hols(1:length(koff_vals),:);
cols = winter(length(koff_vals));
pols = pink(length(koff_vals))/2;
gols = summer(length(koff_vals))/1.5;
vols = parula(4)./2;

%% Dose Response Curves

clear all
f = figure('Color', 'white','Units', 'centimeters','Position',[10 10 37 37],...
    'PaperUnits','centimeters','PaperPosition', [1 2 37 37], 'Renderer','painters');

% Position of axes
left = [0.1, 0.55, 0.1, 0.55];
bottom = [0.55, 0.55, 0.1, 0.1];
width = 0.40;
height = 0.38;

for i=1:4
    ax(i) = axes(f, 'Units','normalized','Position',[left(i) bottom(i) width height]);
end
hold(ax(:),'on')

a1_pos = get(ax(1), 'Position');
a2_pos = get(ax(2), 'Position');
a3_pos = get(ax(3), 'Position');
a4_pos = get(ax(4), 'Position');

% New Simulations
koff_vals = [0.0001, 0.001, 0.01, 0.1];
kd_vals = koff_vals ./ 0.1;
rho_vals = linspace(-4,3.5,16);

kd_lbls = arrayfun(@num2str, kd_vals, 'UniformOutput',0);
xlbl = 10.^(-3:3);
xlabels = arrayfun(@num2str, xlbl, 'UniformOutput',0);

% Plot colors
hols = autumn(length(koff_vals)+1);
hols = hols(1:length(koff_vals),:);
cols = winter(length(koff_vals));

% Axis scale for DR Curves
set(ax(1:4),'xscale','log','tickdir','out')
hold(ax(1:4),'on')

scat1 = DoseResponsePlots(ax(1), kd_vals, rho_vals, 1, 5, 'bound', cols, false, true,'TauLeap');
scat2 = DoseResponsePlots(ax(2), kd_vals, rho_vals, 20, 5, 'bound', hols, false, false,'TauLeap');
scat3 = DoseResponsePlots(ax(3), kd_vals, rho_vals, 1, 5, 'phos', cols, true, true,'TauLeap');
scat4 = DoseResponsePlots(ax(4), kd_vals, rho_vals, 20, 5, 'phos', hols, true, false,'TauLeap');

set(ax(1), 'Position', a1_pos);
set(ax(2), 'Position', a2_pos);
set(ax(3), 'Position', a3_pos);
set(ax(4), 'Position', a4_pos);

ylim(ax(1:4), [0 300]);
xticks(ax(1:4),xlbl);
xticklabels(ax(1:2), [])
xticklabels(ax(3:4), xlabels)
yticklabels(ax([2,4]), [])
set(ax(1:4), 'box', 'off')

% Legend Format
l1 = legend(ax(3), arrayfun(@num2str,kd_vals,'UniformOutput',0),'Units','normalized',...
        'Position',[0.20 0.44 0.1 0.01],'NumColumns',4, 'FontSize', 14);
title(l1, '$K_D$', 'Interpreter','Latex');

l2 = legend(ax(4), arrayfun(@num2str,kd_vals,'UniformOutput',0),'Units','normalized',...
        'Position',[0.65 0.44 0.1 0.01], 'NumColumns',4, 'FontSize', 14);
title(l2, '$K_D$', 'Interpreter','Latex');

title(ax(2), '20 TPC Clusters','Color','r')
title(ax(1), 'Uniform','Color','k')

set(findall(f,'-property','FontSize'),'FontSize',16)
set(findall(f,'-property','Interpreter'),'Interpreter','Latex')
figure(f);

% Save
set(f, 'Units', 'inches')
pos = get(f, 'Position');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
%print(f, '../Figures/Manuscript/Final-Figures/PhosTCR','-dpdf','-r0');

%%
% Axis scales for Coop Plot

%xlim(ax(5), [0 2]);
%ylim(ax(5), [0.001 100]);
%xticks(ax(5),[0.001, 0.01, 0.1, 1]);

%[alpha1, line1] = CooperativityPlots(ax(5), kd_vals, rho_vals, tcrs_per_cluster, valence, bound_type{btype}, hols);
%[alpha2, line2] = CooperativityPlots(ax(5), kd_vals, rho_vals, tcrs_per_cluster2, valence, bound_type{btype}, cols);
%[alpha3, line3] = CooperativityPlots(ax(5), kd_vals, rho_vals, tcrs_per_cluster3, valence, bound_type{btype}, gols);
%[alpha4, line4] = CooperativityPlots(ax(5), kd_vals, rho_vals, tcrs_per_cluster4, valence, bound_type{btype}, pols);

koff_vals = [0.0001, 0.0003, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3];
kd_vals = koff_vals ./ 0.1;

figure()

line1 = Ec50_Emax(gca, kd_vals, rho_vals, tcrs_per_cluster, valence, bound_type{btype}, hols);
line2 = Ec50_Emax(gca, kd_vals, rho_vals, tcrs_per_cluster2, valence, bound_type{btype}, cols);
line3 = Ec50_Emax(gca, kd_vals, rho_vals, tcrs_per_cluster3, valence, bound_type{btype}, gols);
line4 = Ec50_Emax(gca, kd_vals, rho_vals, tcrs_per_cluster4, valence, bound_type{btype}, pols);

set(gca,'yscale','log','tickdir','out')

% Legend
l1 = legend(ax(1), scat1, kd_lbls,'Position',a1_pos.*[1 1 0.1 0.1]+[0.1 0.325 0 0],'Units','normalized','NumColumns',4);
l2 = legend(ax(2), scat2, kd_lbls,'Position',a2_pos.*[1 1 0.1 0.1]+[0.1 0.325 0 0],'Units','normalized','NumColumns',4);
l3 = legend(ax(3), scat3, kd_lbls,'Position',a3_pos.*[1 1 0.1 0.1]+[0.1 0.325 0 0],'Units','normalized','NumColumns',4);
l4 = legend(ax(4), scat4, kd_lbls,'Position',a4_pos.*[1 1 0.1 0.1]+[0.1 0.325 0 0],'Units','normalized','NumColumns',4);

%label1 = ['Uni: ',num2str(round(alpha1,2),'%0.2f')];
%label2 = [num2str(tcrs_per_cluster2),'TPC: ',num2str(round(alpha2,2),'%0.2f')];
%label3 = [num2str(tcrs_per_cluster3),'TPC: ',num2str(round(alpha3,2),'%0.2f')];
%label4 = [num2str(tcrs_per_cluster4),'TPC: ',num2str(round(alpha4,2),'%0.2f')];

label1 = ['Uniform'];
label2 = [num2str(tcrs_per_cluster2),'TPC'];
label3 = [num2str(tcrs_per_cluster3),'TPC'];
label4 = [num2str(tcrs_per_cluster4),'TPC'];

%l5 = legend(ax(5), [line1 line2 line3 line4], {label1, label2, label3, label4}, 'location', 'northwest');

%l1 = legend(ax(1), [scat11 scat21 scat31 scat41], {'1','2','5','10'},'location','northwest','NumColumns',2);
%l2 = legend(ax(2), [scat12 scat22 scat32 scat42], {'1','2','5','10'},'location','northwest','NumColumns',2);

%l3 = legend(ax(3), [line1 line2], {num2str(alpha1), num2str(alpha2)}, 'location', 'northwest');

title(l1, ['Uniform $K_D$']);
title(l2, [num2str(tcrs_per_cluster2),'TPC $K_D$']);
title(l3, [num2str(tcrs_per_cluster3),'TPC $K_D$']);
title(l4, [num2str(tcrs_per_cluster4),'TPC $K_D$']);
title(l5, 'Slope');

%title(l3, 'Slope $\alpha$');

set(findall(f,'-property','FontSize'),'FontSize',14)

% Subplot Letters

for i=1:4
    text(ax(i),0,1.1,charlbls{i},'Units','normalize','FontSize',32);
end

text(ax(5),0,1.05,charlbls{5},'Units','normalize','FontSize',32);

set(findall(f,'-property','Interpreter'),'Interpreter','Latex')
%% Single Dose Response 

bound_type = {'phos', 'bound', 'np', 'covered','fracphos', 'xPMHC', 'rescaled'};
btype = 6;                          % 1 for phos; 2 for bound; 3 for NPs
tcrs_per_cluster = '20'; % ['20x10'];
valence = 10;
kp = [0.1,2];

% Kd and NP rho
%koff_vals = [0.0001, 0.0003, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3];
koff_vals = [0.0001, 0.001, 0.01, 0.1];
kd_vals = koff_vals ./ 0.1;
rho_vals = linspace(-4,3.5,16);
rho_rectangles = linspace(-4.25,3.75,17);
rho_rectangles = 10.^rho_rectangles;

f = figure('Color','white','Position',[10 100 1500 800]);
ax(1) = subplot(221); ax(2) = subplot(222); ax(3) = subplot(223); ax(4) = subplot(224);
a1_pos = get(ax(1), 'Position');
a2_pos = get(ax(2), 'Position');

hols = autumn(length(koff_vals)+1);
hols = hols(1:length(koff_vals),:);
cols = gray(length(rho_vals));

kd_lbls = arrayfun(@num2str, kd_vals, 'UniformOutput',0);
xlbl = 10.^(-3:3);
xlabels = arrayfun(@num2str, xlbl, 'UniformOutput',0);

% Axis scale for DR Curves
set(ax(1),'xscale','log','tickdir','out')
ylim(ax(1), [0 1])

scat1 = DoseResponsePlots(ax(1), kd_vals, kp, rho_vals, tcrs_per_cluster, valence, bound_type{btype}, hols, true, true,'TauLeap');
patch(ax(1), [rho_rectangles(2) rho_rectangles(2) rho_rectangles(3) rho_rectangles(3)],...
    [-10 300 300 -10], cols(1,:), 'FaceAlpha', 0.2, 'EdgeColor','none');

patch(ax(1), [rho_rectangles(7) rho_rectangles(7) rho_rectangles(8) rho_rectangles(8)],...
    [-10 300 300 -10], cols(7,:), 'FaceAlpha', 0.2, 'EdgeColor','none');

patch(ax(1), [rho_rectangles(15) rho_rectangles(15) rho_rectangles(16) rho_rectangles(16)],...
    [-10 300 300 -10], cols(15,:), 'FaceAlpha', 0.2, 'EdgeColor','none');

xticks(ax(1),xlbl);
xticklabels(ax(1), xlabels);
set(ax(1:4), 'box', 'off')
grid(ax(1:4), 'on');

% Axis scales for Coop Plot
%{
set(ax(2),'xscale','log','yscale','log','tickdir','out')
xlim(ax(2), [0.0005 20]);
ylim(ax(2), [0.001 100]);
xticks(ax(2),[0.001, 0.01, 0.1, 1, 10]);

%[alpha1, line1] = CooperativityPlots(ax(2), kd_vals, rho_vals, tcrs_per_cluster, valence, bound_type{btype}, hols);
%}

NoiseScatter(ax(2), kd_vals, rho_vals, tcrs_per_cluster, valence, bound_type{btype}, hols);

%{
for i=1:length(koff_vals)
    hold(ax(3),'on')
    Plot_dprime(ax(3), [koff_vals(i), koff_vals(i)], rho_vals, [valence, valence], tcrs_per_cluster, bound_type{btype}, hols(i,:), true)
end

rho_vals = rho_vals([2, 7, 11, 15]);

for i=1:length(rho_vals)
    hold(ax(4),'on')
    Plot_dprime(ax(4), koff_vals, [rho_vals(i), rho_vals(i)], [valence, valence], tcrs_per_cluster, bound_type{btype}, cols(i,:), false)
end

%Plot_dprime(ax(4), koff_vals, [rho_vals(8), rho_vals(8)], [valence, valence], 3, bound_type{btype}, cols, false)
%}

SpecificityPlots([ax(3) ax(4)], kd_vals, rho_vals, tcrs_per_cluster, valence, bound_type{btype}, hols)

set(ax(3:4),'xscale','log','tickdir','out','box','off');
grid(ax(3:4), 'on');

set(ax(1), 'Position', a1_pos);
set(ax(2), 'Position', a2_pos);

%l1 = legend(ax(2), line1, ['Slope = ', num2str(alpha1)], 'Location', 'southeast');
%l2 = legend(ax(4), arrayfun(@num2str, rho_vals, 'UniformOutput', false), 'Location', 'northeast', 'NumColumns', 4);
%title(l1, ['$K_D$'])

rho_lbls = linspace(-4,3,8);
rho_lbls = arrayfun(@num2str, rho_lbls, 'UniformOutput',0);
rho_lbls = append('$10^{',rho_lbls);
rho_lbls = append(rho_lbls,'}$');

colormap(ax(3), autumn)
colormap(ax(4), autumn)
c = colorbar(ax(3:4), 'Location', 'eastoutside', 'Ticks', linspace(0,1,length(kd_vals)),...
    'TickLabels', {'$10^{-3}$', '', '$10^{-2}$', '', '$10^{-1}$', '', '$10^0$', '', '$10^1$'}, 'TickLabelInterpreter','latex');

title(c,'$K_D$', 'Interpreter', 'latex');


set(findall(f,'-property','FontSize'),'FontSize',14)

set(findall(f,'-property','Interpreter'),'Interpreter','Latex')

%% Dose Response to Valence

f = figure('Color', 'white','Position',[10 100 1500 1000]);
hols = autumn(7);
cols = winter(7);
pols = pink(7)./1.5;
gols = summer(7);
vols = parula(4)./2;

charlbls = {'A','B','C','D','E'};

tcrs_per_cluster = 3;
tcrs_per_cluster2 = 3;
tcrs_per_cluster3 = 5;
tcrs_per_cluster4 = 20;
valence = [1,2,3,4,5,10];
bound_type = {'phos', 'bound', 'np', 'covered','fracphos','normBound'};
btype = 6;

% New Simulations
koff_vals = [0.0001, 0.001, 0.01, 0.1];
kd_vals = koff_vals ./ 0.1;
rho_vals = linspace(-4,3.5,16);

ax(1) = subplot(2,2,1); ax(2) = subplot(2,2,2); ax(3)=subplot(2,2,3); ax(4) = subplot(2,2,4);

set(ax(1:4),'xscale','log','tickdir','out')
ylim(ax(1:4), [0 300]);

a1_pos = get(ax(1), 'Position');
a2_pos = get(ax(2), 'Position');
a3_pos = get(ax(3), 'Position');
a4_pos = get(ax(4), 'Position');

for i=1:length(valence)
    clus1(i) = DoseResponsePlots(ax(1), kd_vals(1), rho_vals, tcrs_per_cluster,...
        valence(i), bound_type{btype}, pols(8-i,:), false, true, 'TauLeap');
    clus2(i) = DoseResponsePlots(ax(2), kd_vals(2), rho_vals, tcrs_per_cluster,...
        valence(i), bound_type{btype}, pols(8-i,:), false, false, 'TauLeap');
    clus3(i) = DoseResponsePlots(ax(3), kd_vals(3), rho_vals, tcrs_per_cluster,...
        valence(i), bound_type{btype}, pols(8-i,:), true, true, 'TauLeap');
    clus4(i) = DoseResponsePlots(ax(4), kd_vals(4), rho_vals, tcrs_per_cluster,...
        valence(i), bound_type{btype}, pols(8-i,:), true, false, 'TauLeap');

end
%uni = DoseResponsePlots(ax(1), kd_vals, rho_vals, tcrs_per_cluster, valence(1), bound_type{2}, hols(1,:), true, true);

set(ax(1), 'Position', a1_pos);
set(ax(2), 'Position', a2_pos);
set(ax(3), 'Position', a3_pos);
set(ax(4), 'Position', a4_pos);

% Axis scale for DR Curves
xlbl = 10.^(-3:3);
xlabels = arrayfun(@num2str, xlbl, 'UniformOutput',0);
xticks(ax(1:4),xlbl);
xticklabels(ax(1:4), xlabels)
set(ax(1:4), 'box', 'off')

%l1 = legend(ax(1), uni, '1-10','Location','northwest','NumColumns',3);
l1 = legend(ax(1), clus1, arrayfun(@num2str, valence,'UniformOutput',0),'location','northwest','NumColumns',3);
l2 = legend(ax(2), clus2, arrayfun(@num2str, valence,'UniformOutput',0),'location','northwest','NumColumns',3);
l3 = legend(ax(3), clus3, arrayfun(@num2str, valence,'UniformOutput',0),'location','northwest','NumColumns',3);
l4 = legend(ax(4), clus4, arrayfun(@num2str, valence,'UniformOutput',0),'location','northwest','NumColumns',3);

title(l1,['Valence || $K_D$ ', num2str(kd_vals(1))])
title(l2,['Valence || $K_D$ ', num2str(kd_vals(2))])
title(l3,['Valence || $K_D$ ', num2str(kd_vals(3))])
title(l4,['Valence || $K_D$ ', num2str(kd_vals(4))])

set(findall(f,'-property','FontSize'),'FontSize',14)

for i=1:4
    text(ax(i),0,1.1,charlbls{i},'Units','normalize','FontSize',32);
end

set(findall(f,'-property','Interpreter'),'Interpreter','Latex')

%% Surface Comparison

clear all
f1 = figure('Color', 'white','Units', 'centimeters','Position',[1 1 45 20],...
    'PaperUnits','centimeters','PaperPosition', [1 2 45 20], 'Renderer','painters','visible','off');

f2 = figure('Color', 'white','Units', 'centimeters','Position',[1 1 45 20],...
    'PaperUnits','centimeters','PaperPosition', [1 2 45 20], 'Renderer','painters','visible','off');

% Position of axes
left = [0.1, 0.37, 0.63, 0.1, 0.37, 0.63];
bottom = [0.6, 0.6, 0.6, 0.1, 0.1, 0.1];
width = 0.23;
height = 0.35;

for i=1:6
    ax(i) = axes(f1, 'Units','normalized','Position',[left(i) bottom(i) width height]);
end

hold(ax(:),'on')

% New Simulations
koff_vals = [0.0001, 0.001, 0.01];
kd_vals = koff_vals ./ 0.1;
rho_vals = linspace(-4,3.5,16);
bound_type = {'phos', 'bound', 'np', 'covered','fracphos', 'normBound', 'normPhos','rescaled'};

valence = [5];
tpc = [1, 3, 10];
kp = [0.01,2];
np_radius = 20;

% Color and Linestyle
col = gray(length(valence)+2);
hol = winter(length(valence)+2);
sol = autumn(length(valence)+2);

linestyle = {'-','--'};
markerstyle = {'s','d'};
colorstyle = {hol, hol, sol};

% Kd Titles
for i=1:3
     title(ax(i),['$K_{D} = ',num2str(kd_vals(i)),'$']);
end

for i=1:length(kd_vals)
    for j=1:length(valence)
        for k=1:length(tpc)
            utop{j,k} = DoseResponsePlots(ax(i), kd_vals(i), kp, rho_vals, tpc(k),...
            valence(j), bound_type{6}, colorstyle{k}(j,:), false, false, np_radius);
    
            ubot{j,k} = DoseResponsePlots(ax(i+3), kd_vals(i), kp, rho_vals, tpc(k),...
            valence(j), bound_type{7}, colorstyle{k}(j,:), false, false, np_radius);
        end
    end
end

%koff_vals = [0.0001, 0.0003, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3];
koff_vals = [0.0001, 0.001, 0.01, 0.1];
kd_vals = koff_vals ./ 0.1;

% Cooperativity axes positions
l2 = [0.1, 0.59];
b2 = [0.2, 0.2];
w2 = 0.40;
h2 = 0.75;

for i=7:8
    ax(i) = axes(f2, 'Units','normalized','Position',[l2(mod(i,6)) b2(mod(i,6)) w2 h2]);
end
ax(9) = axes(f2, 'Units','normalized','Position',[0.30 0.35 0.15 0.18]);

hold(ax(:),'on')

for j=1:length(valence)
    for k=1:length(tpc)
        [a1, ln1] = CooperativityPlots(ax(7), kd_vals, kp, rho_vals, tpc(k), valence(j), np_radius,...
            bound_type{6}, colorstyle{k}(j,:),true,true,1);
    
        [a4, ln4] = CooperativityPlots(ax(8), kd_vals, kp, rho_vals, tpc(k), valence(j), np_radius,...
            bound_type{6}, colorstyle{k}(j,:),true,true,2);
    
        hold(ax(9),'on')
        scatter(ax(9), tpc(k), a1, 72, colorstyle{k}(j,:),'x');
    end
end

% Text fontsize
set(findall([f1 f2],'-property','FontSize'),'FontSize',20)
set(findall([f1 f2],'-property','Interpreter'),'Interpreter','Latex')

% Format First Row
set(ax(1:8),'xscale','log','tickdir','out','box','off')
ylim(ax(1:3), [0 300]);
xlim(ax(1:6), [0.0001 15000]);
xlabel(ax(2), 'pMHC Concentration (log$_{10}$ a.u.)');
ylabel(ax(1), 'Bound TCRs')

% Format Second Row
ylim(ax(4:6), [0 300]);
yticks(ax(4:6), [0 100 200 300]);
%xlim(ax(4:6), [0.0001 2000]);
xlabel(ax(5), 'pMHC Concentration (log$_{10}$ a.u.)');
ylabel(ax(4), ' Phos. TCRs');

% Format Third Row
set(ax(7),'yscale','log');

xlim(ax(7:8), [0.0007 2]);
ylim(ax(7), [0.0001 100]);
ylim(ax(8), [0 300]);

xlbl = 10.^(-3:1:0);
xlabels = {'-3', '-2', '-1', '0'};
xticks(ax(7:8), xlbl);
xticklabels(ax(7:8), xlabels);

ylbl = 10.^(-4:2:2);
ylabels = {'-4','-2','0','2'};
yticks(ax(7), ylbl);
yticklabels(ax(7), ylabels);

xlabel(ax(7:8), '$K_D$ (log$_{10}$ a.u.)', 'FontSize', 48)
ylabel(ax(7), 'EC50 (log$_{10}$)', 'FontSize', 48)
ylabel(ax(8), 'EMax', 'FontSize', 48)

% Format Inset Plot
grid(ax(9),'on');
set(ax(9), 'box','on','FontSize',16)
xlim(ax(9),[0 21]);
ylim(ax(9),[0.4 1.5])
xlabel(ax(9), 'TPC','FontSize',24)
ylabel(ax(9), 'Slope','FontSize',24)

% Axis Ticks and Format
set(ax(1:8),'GridLineWidth',0.5,'Linewidth',2)
xlbl = 10.^(-3:2:3);
xlabels = {'-3', '-1', '1', '3'};
xticks(ax(1:6),xlbl);
xticklabels(ax(1:6), xlabels)
xsecondarylabel([ax(2) ax(5)], '')
set(ax(1:6),'TickLabelInterpreter', 'Latex')
yticklabels(ax(2:3), '')
yticklabels(ax(5:6), '')

% Legend Format
l1 = legend(ax(8), [utop{1,:}], arrayfun(@num2str,tpc,'UniformOutput',0),'Units','normalized',...
        'Position',[0.88 0.69 0.1 0.15], 'FontSize', 20);
title(l1, 'TPC');

% Correct axes position
for i=1:6
    set(ax(i),'Units','normalized','Position',[left(i) bottom(i) width height]);
end

%{
charlbls = {'A','B','C','D','E','F','G','H'};
for i=1:8
    text(ax(i),0.01,0.9,charlbls{i},'Units','normalize','FontSize',24);
end
%}


%set(findall([f1 f2],'-property','FontSize'),'FontSize',48)
set(findall([f1 f2],'-property','LineWidth'),'LineWidth',2)

set(findall(f1,'-property','Interpreter'),'Interpreter','Latex')
set(findall(f2,'-property','Interpreter'),'Interpreter','Latex')


figure(f1); figure(f2);

% Save
set(f1, 'Units', 'inches')
pos = get(f1, 'Position');
set(f1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
%print(f, '../Figures/Manuscript/Final-Figures/KPR1','-dpdf','-r0');


%% Valence Comparison: 1TPC Uniform

clear all
f = figure('Color', 'white','Units', 'centimeters','Position',[10 10 45 20],...
    'PaperUnits','centimeters','PaperPosition', [1 2 45 20], 'Renderer','painters','visible','off');
f2 = figure('Color', 'white','Units', 'centimeters','Position',[10 10 45 20],...
    'PaperUnits','centimeters','PaperPosition', [1 2 45 20], 'Renderer','painters','visible','off');

% Position of axes
left = [0.1, 0.37, 0.63, 0.1, 0.37, 0.63];
bottom = [0.81, 0.81, 0.81, 0.57, 0.57, 0.57];
width = 0.23;
height = 0.15;

for i=1:6
    ax(i) = axes(f, 'Units','normalized','Position',[left(i) bottom(i) width height]);
end
hold(ax(:),'on')

% New Simulations
koff_vals = [0.0001, 0.001, 0.01];
kd_vals = koff_vals ./ 0.1;
rho_vals = linspace(-4,3.5,16);
bound_type = {'phos', 'bound', 'np', 'covered','fracphos', 'normBound', 'normPhos','rescaled'};

valence = [1,2,3,4,5,10];
tpc = 1;
kp = [0.1,2];

% Colormaps
col = gray(length(valence)+2);
hol = winter(length(kd_vals)+2);
sol = autumn(length(valence)+2);

% Kd Titles
for i=1:3
     title(ax(i),['$K_{D} = ',num2str(kd_vals(i)),'$']);
end

for i=1:length(kd_vals)
    for j=1:length(valence)
        %hold(ax(i),'on');
        utop{j} = DoseResponsePlots(ax(i), kd_vals(i), kp, rho_vals, tpc,...
        valence(j), bound_type{6}, col(j,:), false, false, 'TauLeap');
    
        ubot{j} = DoseResponsePlots(ax(i+3), kd_vals(i), kp, rho_vals, tpc,...
        valence(j), bound_type{3}, col(j,:), false, false, 'TauLeap');
    end
end

%koff_vals = [0.0001, 0.0003, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3];
koff_vals = [0.0001, 0.001, 0.01, 0.1];
kd_vals = koff_vals ./ 0.1;

% Cooperativity axes positions
l2 = [0.1, 0.59];
b2 = [0.2, 0.2];
w2 = 0.40;
h2 = 0.75;

for i=7:8
    ax(i) = axes(f2, 'Units','normalized','Position',[l2(mod(i,6)) b2(mod(i,6)) w2 h2]);
end
ax(9) = axes(f2, 'Units','normalized','Position',[0.30 0.3 0.15 0.18]);

hold(ax(:),'on')

for j=1:length(valence)
    [a1, ln1] = CooperativityPlots(ax(7), kd_vals, kp, rho_vals, tpc, valence(j),...
        bound_type{6}, col(j,:),true,true,1, 'TauLeap');
    
    [a4, ln4] = CooperativityPlots(ax(8), kd_vals, kp, rho_vals, tpc, valence(j),...
        bound_type{6}, col(j,:),true,true,2, 'TauLeap');
    
    hold(ax(9),'on')
    scatter(ax(9), valence(j), a1, 72, col(j,:),'x');
end

% Text fontsize
set(findall([f f2],'-property','FontSize'),'FontSize',20)
set(findall([f f2],'-property','Interpreter'),'Interpreter','Latex')

% Format First Row
set(ax(1:8),'xscale','log','tickdir','out','box','off')
ylim(ax(1:3), [0 300]);
xlim(ax(1:3), [0.0001 15000]);
xlabel(ax(2), 'pMHC Concentration (log$_{10}$ a.u.)');
ylabel(ax(1), 'Bound TCRs')

% Format Second Row
ylim(ax(4:6), [0 200]);
yticks(ax(4:6), [0 100 200]);
xlim(ax(4:6), [0.0001 2000]);
xlabel(ax(5), 'NP Concentration (log$_{10}$ a.u.)');
ylabel(ax(4), ' Bound NPs');

% Format Third Row
set(ax(7),'yscale','log');
xlim(ax(7:8), [0.0007 2]);
ylim(ax(7), [0.001 100]);
ylim(ax(8), [0 250]);

xlbl = 10.^(-3:1:0);
xlabels = {'-3', '-2', '-1', '0'};
xticks(ax(7:8), xlbl);
xticklabels(ax(7:8), xlabels);

ylbl = 10.^(-4:2:2);
ylabels = {'-4','-2','0','2'};
yticks(ax(7), ylbl);
yticklabels(ax(7), ylabels);

xlabel(ax(7:8), '$K_D$ (log$_{10}$ a.u.)', 'FontSize', 48)
ylabel(ax(7), 'EC50 (log$_{10}$)', 'FontSize', 48)
ylabel(ax(8), 'EMax', 'FontSize', 48)

% Format Inset Plot
grid(ax(9),'on');
set(ax(9), 'box','on','FontSize',16)
xlim(ax(9),[0 11]);
ylim(ax(9),[0.3 0.7])
xlabel(ax(9), 'Valence','FontSize',24)
ylabel(ax(9), 'Slope','FontSize',24)

% Axis Ticks and Format
set(ax(1:8),'GridLineWidth',0.5,'Linewidth',2)
xlbl = 10.^(-3:2:3);
xlabels = {'-3', '-1', '1', '3'};
xticks(ax(1:6),xlbl);
xticklabels(ax(1:6), xlabels)
xsecondarylabel([ax(2) ax(5)], '')
set(ax(1:6),'TickLabelInterpreter', 'Latex')
yticklabels(ax(2:3), '')
yticklabels(ax(5:6), '')

% Legend Format

l1 = legend(ax(8), [utop{:}], arrayfun(@num2str,valence,'UniformOutput',0),'Units','normalized',...
        'Position',[0.88 0.69 0.1 0.15], 'FontSize', 20);
title(l1, 'NP Valence');

% Correct axes position
for i=1:6
    set(ax(i),'Units','normalized','Position',[left(i) bottom(i) width height]);
end

%{
charlbls = {'A','B','C','D','E','F','G','H'};
for i=1:8
    text(ax(i),0.01,0.9,charlbls{i},'Units','normalize','FontSize',24);
end
%}

%set(findall(f,'-property','FontSize'),'FontSize',20)
set(findall([f f2],'-property','LineWidth'),'LineWidth',2)
set(findall([f f2],'-property','Interpreter'),'Interpreter','Latex')

figure(f); figure(f2);

% Save
set(f, 'Units', 'inches')
pos = get(f, 'Position');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
%print(f, '../Figures/Manuscript/Final-Figures/F2','-dpdf','-r0');

%% Valence Comparison: 20 TPC Cluster

clear all
f = figure('Color', 'white','Units', 'centimeters','Position',[10 10 45 20],...
    'PaperUnits','centimeters','PaperPosition', [1 2 45 20], 'Renderer','painters','visible','off');
f2 = figure('Color', 'white','Units', 'centimeters','Position',[10 10 45 20],...
    'PaperUnits','centimeters','PaperPosition', [1 2 45 20], 'Renderer','painters','visible','off');

% Position of axes
left = [0.1, 0.37, 0.63, 0.1, 0.37, 0.63];
bottom = [0.81, 0.81, 0.81, 0.57, 0.57, 0.57];
width = 0.23;
height = 0.15;

for i=1:6
    ax(i) = axes(f, 'Units','normalized','Position',[left(i) bottom(i) width height]);
end
hold(ax(:),'on')

% New Simulations
koff_vals = [0.0001, 0.001, 0.01];
kd_vals = koff_vals ./ 0.1;
rho_vals = linspace(-4,3.5,16);
bound_type = {'phos', 'bound', 'np', 'covered','fracphos', 'normBound', 'normPhos','rescaled'};

valence = [1,2,3,4,5,10];
tpc = 20;
kp=[0.1,2];

% Colormaps
col = hot(length(valence)*2);
hol = winter(length(kd_vals)+2);
sol = autumn(length(valence)+2);

% Kd Titles
for i=1:3
     title(ax(i),['$K_{D} = ',num2str(kd_vals(i)),'$']);
end

for i=1:length(kd_vals)
    for j=1:length(valence)
        %hold(ax(i),'on');
        utop{j} = DoseResponsePlots(ax(i), kd_vals(i), kp, rho_vals, tpc,...
        valence(j), bound_type{6}, col(j,:), false, false, 'TauLeap');
    
        ubot{j} = DoseResponsePlots(ax(i+3), kd_vals(i), kp, rho_vals, tpc,...
        valence(j), bound_type{3}, col(j,:), false, false, 'TauLeap');
    end
end

%koff_vals = [0.0001, 0.0003, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3];
koff_vals = [0.0001, 0.001, 0.01, 0.1];
kd_vals = koff_vals ./ 0.1;

% Cooperativity axes positions
l2 = [0.1, 0.59];
b2 = [0.2, 0.2];
w2 = 0.40;
h2 = 0.75;

for i=7:8
    ax(i) = axes(f2, 'Units','normalized','Position',[l2(mod(i,6)) b2(mod(i,6)) w2 h2]);
end
ax(9) = axes(f2, 'Units','normalized','Position',[0.35 0.29 0.15 0.18]);

hold(ax(:),'on')

for j=1:length(valence)
    [a1, ln1] = CooperativityPlots(ax(7), kd_vals, kp, rho_vals, tpc, valence(j),...
        bound_type{6}, col(j,:),true,true,1, 'TauLeap');
    
    [a4, ln4] = CooperativityPlots(ax(8), kd_vals, kp, rho_vals, tpc, valence(j),...
        bound_type{6}, col(j,:),true,true,2, 'TauLeap');
    
    hold(ax(9),'on')
    scatter(ax(9), valence(j), a1, 72, col(j,:),'x');
end

% Text fontsize
set(findall([f f2],'-property','FontSize'),'FontSize',20)
set(findall([f f2],'-property','Interpreter'),'Interpreter','Latex')

% Format First Row
set(ax(1:8),'xscale','log','tickdir','out','box','off')
ylim(ax(1:3), [0 300]);
xlim(ax(1:3), [0.0001 15000]);
xlabel(ax(2), 'pMHC Concentration (log$_{10}$ a.u.)');
ylabel(ax(1), 'Bound TCRs')

% Format Second Row
ylim(ax(4:6), [0 200]);
yticks(ax(4:6), [0 100 200]);
xlim(ax(4:6), [0.0001 2000]);
xlabel(ax(5), 'NP Concentration (log$_{10}$ a.u.)');
ylabel(ax(4), ' Bound NPs');

% Format Third Row
set(ax(7),'yscale','log');
xlim(ax(7:8), [0.0007 2]);
ylim(ax(7), [0.001 100]); %[0.001 100]
ylim(ax(8), [0 250]);

xlbl = 10.^(-2:1:0);
xlabels = {'-2', '-1', '0'};
xticks(ax(7:8), xlbl);
xticklabels(ax(7:8), xlabels);

ylbl = 10.^(-4:2:2);
ylabels = {'-4','-2','0','2'};
yticks(ax(7), ylbl);
yticklabels(ax(7), ylabels);

xlabel(ax(7:8), '$K_D$ (log$_{10}$ a.u.)', 'FontSize', 48)
ylabel(ax(7), 'EC50 (log$_{10}$)', 'FontSize', 48)
ylabel(ax(8), 'EMax', 'FontSize', 48)

% Format Inset Plot
grid(ax(9),'on');
set(ax(9), 'box','on','FontSize',16)
xlim(ax(9),[0 11]);
ylim(ax(9),[0.3 1.5])
xlabel(ax(9), 'Valence','FontSize',24)
ylabel(ax(9), 'Slope','FontSize',24)

% Axis Ticks and Format
set(ax(1:8),'GridLineWidth',0.5,'Linewidth',2)
xlbl = 10.^(-3:2:3);
xlabels = {'-3', '-1', '1', '3'};
xticks(ax(1:6),xlbl);
xticklabels(ax(1:6), xlabels)
xsecondarylabel([ax(2) ax(5)], '')
set(ax(1:6),'TickLabelInterpreter', 'Latex')
yticklabels(ax(2:3), '')
yticklabels(ax(5:6), '')

% Legend Format
l1 = legend(ax(8), [utop{:}], arrayfun(@num2str,valence,'UniformOutput',0),'Units','normalized',...
        'Position',[0.88 0.69 0.1 0.15], 'FontSize', 20);
title(l1, 'NP Valence');

% Correct axes position
for i=1:6
    set(ax(i),'Units','normalized','Position',[left(i) bottom(i) width height]);
end

%{
charlbls = {'A','B','C','D','E','F','G','H'};
for i=1:8
    text(ax(i),0.01,0.9,charlbls{i},'Units','normalize','FontSize',24);
end
%}

%set(findall(f,'-property','FontSize'),'FontSize',14)
set(findall([f f2],'-property','LineWidth'),'LineWidth',2)

set(findall([f f2],'-property','Interpreter'),'Interpreter','Latex')
figure(f);figure(f2);

% Save
set(f, 'Units', 'inches')
pos = get(f, 'Position');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
%print(f, '../Figures/Manuscript/Final-Figures/F3','-dpdf','-r0');

%% r14 to r20 Comparison

clear all
f = figure('Color', 'white','Units', 'centimeters','Position',[10 10 37 37],...
    'PaperUnits','centimeters','PaperPosition', [1 2 37 37], 'Renderer','painters','visible','off');

% Position of axes
left = [0.1, 0.37, 0.63, 0.1, 0.37, 0.63];
bottom = [0.81, 0.81, 0.81, 0.57, 0.57, 0.57];
width = 0.23;
height = 0.15;

for i=1:6
    ax(i) = axes(f, 'Units','normalized','Position',[left(i) bottom(i) width height]);
end
hold(ax(:),'on')

% New Simulations
koff_vals = [0.0001, 0.001, 0.01];
kd_vals = koff_vals ./ 0.1;
rho_vals = linspace(-4,3.5,16);
bound_type = {'phos', 'bound', 'np', 'covered','fracphos', 'normBound', 'normPhos','rescaled'};

% Colormaps
col = cool(6);

tpc = 20;
kp=[0.1,2];

% Kd Titles
for i=1:3
     title(ax(i),['$K_{D} = ',num2str(kd_vals(i)),'$']);
end

for i=1:length(kd_vals)
    utop{1} = DoseResponsePlots(ax(i), kd_vals(i), kp, rho_vals, tpc,...
    2, bound_type{6}, col(1,:), false, false, 'TauLeap');

    c20top{4} = DoseResponsePlots(ax(i), kd_vals(i), kp, rho_vals, tpc,...
    5, bound_type{6}, col(2,:), false, false, 'TauLeap');

    c20top{1} = DoseResponsePlots(ax(i), kd_vals(i), kp, rho_vals, tpc,...
    10, bound_type{6}, col(3,:), false, false, 'TauLeap');

    c20top{2} = DoseResponsePlots(ax(i), kd_vals(i), kp, rho_vals, tpc,...
    1, bound_type{6}, col(5,:), false, false, 'r14');

    c20top{3} = DoseResponsePlots(ax(i), kd_vals(i), kp, rho_vals, tpc,...
    5, bound_type{6}, col(6,:), false, false, 'r14');

    ubot{1} = DoseResponsePlots(ax(i+3), kd_vals(i), kp, rho_vals, tpc,...
    2, bound_type{3}, col(1,:), false, false, 'TauLeap');

    c20bot{4} = DoseResponsePlots(ax(i+3), kd_vals(i), kp, rho_vals, tpc,...
    5, bound_type{3}, col(2,:), false, false, 'TauLeap');

    c20bot{1} = DoseResponsePlots(ax(i+3), kd_vals(i), kp, rho_vals, tpc,...
    10, bound_type{3}, col(3,:), false, false, 'TauLeap');

    c20bot{2} = DoseResponsePlots(ax(i+3), kd_vals(i), kp, rho_vals, tpc,...
    1, bound_type{3}, col(5,:), false, false, 'r14');

    c20bot{3} = DoseResponsePlots(ax(i+3), kd_vals(i), kp, rho_vals, tpc,...
    5, bound_type{3}, col(6,:), false, false, 'r14');
end

koff_vals = [0.0001, 0.0003, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3];
koff_vals = [0.0001, 0.001, 0.01, 0.1];
kd_vals = koff_vals ./ 0.1;

% Cooperativity axes positions
l2 = [0.1, 0.59];
b2 = [0.1, 0.1];
w2 = 0.40;
h2 = 0.30;

for i=7:8
    ax(i) = axes(f, 'Units','normalized','Position',[l2(mod(i,6)) b2(mod(i,6)) w2 h2]);
end
ax(9) = axes(f, 'Units','normalized','Position',[0.36 0.15 0.12 0.08]);

hold(ax(:),'on')

pmhc_density = [2, 5, 10, 1, 5] ./ [20, 20, 20, 14, 14].^2;

for i=0:1
    [a1{i+1}, ln1] = CooperativityPlots(ax(7+i), kd_vals, kp, rho_vals, tpc, 2,...
        bound_type{6}, col(1,:),true,true,1+i, 'TauLeap');
    
    [a2{i+1}, ln4] = CooperativityPlots(ax(7+i), kd_vals, kp, rho_vals, tpc, 5,...
        bound_type{6}, col(2,:),true,true,1+i, 'TauLeap');

    [a3{i+1}, ln4] = CooperativityPlots(ax(7+i), kd_vals, kp, rho_vals, tpc, 10,...
        bound_type{6}, col(3,:),true,true,1+i, 'TauLeap');

    [a4{i+1}, ln1] = CooperativityPlots(ax(7+i), kd_vals, kp, rho_vals, tpc, 1,...
        bound_type{6}, col(5,:),true,true,1+i, 'r14');

    [a5{i+1}, ln4] = CooperativityPlots(ax(7+i), kd_vals, kp, rho_vals, tpc, 5,...
        bound_type{6}, col(6,:),true,true,1+i, 'r14');
end

hold(ax(9),'on')
scatter(ax(9), pmhc_density, [a1{1} a2{1} a3{1} a4{1} a5{1}], [], [col(1:3,:); col(5:6,:)],'x');

% Text fontsize
set(findall(f,'-property','FontSize'),'FontSize',20)
set(findall(f,'-property','Interpreter'),'Interpreter','Latex')

% Format First Row
set(ax(1:8),'xscale','log','tickdir','out','box','off')
ylim(ax(1:3), [0 300]);
xlim(ax(1:3), [0.0001 15000]);
xlabel(ax(2), 'pMHC Concentration (log$_{10}$ a.u.)');
ylabel(ax(1), 'Bound TCRs')

% Format Second Row
ylim(ax(4:6), [0 200]);
yticks(ax(4:6), [0 100 200]);
xlim(ax(4:6), [0.0001 2000]);
xlabel(ax(5), 'NP Concentration (log$_{10}$ a.u.)');
ylabel(ax(4), ' Bound NPs');

% Format Third Row
set(ax(7),'yscale','log');
xlim(ax(7:8), [0.0007 2]);
ylim(ax(7), [0.00005 100]);
ylim(ax(8), [0 250]);
yticks(ax(7), [0.0001 0.01 1 100]);

xlbl = 10.^(-2:1:0);
xlabels = {'-2', '-1', '0'};
xticks(ax(7:8), xlbl);
xticklabels(ax(7:8), xlabels);

ylbl = 10.^(-4:2:2);
ylabels = {'-4','-2','0','2'};
yticks(ax(7), ylbl);
yticklabels(ax(7), ylabels);

xlabel(ax(7:8), '$K_D$ (log$_{10}$ a.u.)', 'FontSize', 24)
ylabel(ax(7), 'EC50 (log$_{10}$)', 'FontSize', 24)
ylabel(ax(8), 'EMax', 'FontSize', 24)

% Format Inset Plot
grid(ax(9),'on');
set(ax(9), 'box','on','FontSize',10)
xlim(ax(9),[0 0.03]);
ylim(ax(9),[0 1.5])
xlabel(ax(9), 'pMHC Density','FontSize',14)
ylabel(ax(9), 'Slope','FontSize',14)

% Axis Ticks and Format
set(ax(1:8),'GridLineWidth',0.5,'Linewidth',2)
xlbl = 10.^(-3:2:3);
xlabels = {'-3', '-1', '1', '3'};
xticks(ax(1:6),xlbl);
xticklabels(ax(1:6), xlabels)
xsecondarylabel([ax(2) ax(5)], '')
set(ax(1:6),'TickLabelInterpreter', 'Latex')
yticklabels(ax(2:3),'')
yticklabels(ax(5:6),'')

% Legend Format
l1 = legend(ax(8), [utop{1} c20top{4} c20top{1} c20top{2} c20top{3}],...
    {'(20,2)','(20,5)', '(20,10)', '(14,1)','(14,5)'},'Units','normalized',...
        'Position',[0.88 0.69 0.1 0.15], 'FontSize', 20);
title(l1, 'NPs $(r,v)$');

% Correct axes position
for i=1:6
    set(ax(i),'Units','normalized','Position',[left(i) bottom(i) width height]);
end

%set(findall(f,'-property','FontSize'),'FontSize',14)
set(findall(f,'-property','Interpreter'),'Interpreter','Latex')

figure(f);

% Save
set(f, 'Units', 'inches')
pos = get(f, 'Position');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(f, '../Figures/Manuscript/Final-Figures/F4','-dpdf','-r0');

%% Kinetic Proofreading

clear all
f = figure('Color', 'white','Units', 'centimeters','Position',[10 10 37 37],...
    'PaperUnits','centimeters','PaperPosition', [1 2 37 37], 'Renderer','painters','visible','off');

% Position of axes
left = [0.1, 0.37, 0.63, 0.1, 0.37, 0.63];
bottom = [0.81, 0.81, 0.81, 0.57, 0.57, 0.57];
width = 0.23;
height = 0.15;

for i=1:6
    ax(i) = axes(f, 'Units','normalized','Position',[left(i) bottom(i) width height]);
end
hold(ax(:),'on')

% New Simulations
koff_vals = [0.0001, 0.001, 0.01];
kd_vals = koff_vals ./ 0.1;
rho_vals = linspace(-4,3.5,16);
bound_type = {'phos', 'bound', 'np', 'covered','fracphos', 'normBound', 'normPhos','rescaled'};

valence = [5];
tpc = [20];
kp = [0.001, 0.01, 0.1, 1];
n = [2,3,4,5];

% Color and Linestyle
col = gray(length(kp)+2);
hol = winter(length(kp)+2);
sol = autumn(length(kp)+2);

linestyle = {'-','--'};
markerstyle = {'s','d'};
colorstyle = {col, hol, sol};

% Kd Titles
for i=1:3
     title(ax(i),['$K_{D} = ',num2str(kd_vals(i)),'$']);
end

for i=1:length(kd_vals)
    for j=1:length(kp)
        for k=1:length(tpc)
            utop{j,k} = DoseResponsePlots(ax(i), kd_vals(i), [kp(j),n(1)], rho_vals, 1,...
            valence, bound_type{7}, colorstyle{1}(j,:), false, false, 'TauLeap');
    
            ubot{j,k} = DoseResponsePlots(ax(i+3), kd_vals(i), [kp(j),n(1)], rho_vals, 20,...
            valence, bound_type{7}, colorstyle{3}(j,:), false, false, 'TauLeap');
        end
    end
end

%koff_vals = [0.0001, 0.0003, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3];
koff_vals = [0.0001, 0.001, 0.01, 0.1];
kd_vals = koff_vals ./ 0.1;

% Cooperativity axes positions
l2 = [0.1, 0.59];
b2 = [0.1, 0.1];
w2 = 0.40;
h2 = 0.30;

for i=7:8
    ax(i) = axes(f, 'Units','normalized','Position',[l2(mod(i,6)) b2(mod(i,6)) w2 h2]);
end

hold(ax(:),'on')

for j=1:length(kp)
    for k=1:length(tpc)
        [a1, ln1] = CooperativityPlots(ax(7), kd_vals, [kp(j),n(1)], rho_vals, 1, valence,...
            bound_type{7}, colorstyle{1}(j,:),true,true,2, 'TauLeap');
    
        [a4, ln4] = CooperativityPlots(ax(8), kd_vals, [kp(j),n(1)], rho_vals, 20, valence,...
            bound_type{7}, colorstyle{3}(j,:),true,true,2, 'TauLeap');
    
    end
end

% Text fontsize
set(findall(f,'-property','FontSize'),'FontSize',20)
set(findall(f,'-property','Interpreter'),'Interpreter','Latex')

% Format First Row
set(ax(1:8),'xscale','log','tickdir','out','box','off')
ylim(ax(1:3), [0 300]);
xlim(ax(1:6), [0.0001 15000]);
xlabel(ax(2), 'pMHC Concentration (log$_{10}$ a.u.)','FontSize',24);
ylabel(ax(1), 'Phos. TCRs', 'FontSize', 24)

% Format Second Row
ylim(ax(4:6), [0 300]);
yticks(ax(4:6), [0 100 200 300]);
%xlim(ax(4:6), [0.0001 2000]);
xlabel(ax(5), 'pMHC Concentration (log$_{10}$ a.u.)','FontSize',24);
ylabel(ax(4), ' Phos. TCRs','FontSize',24);

% Format Third Row
xlim(ax(7:8), [0.0007 2]);
ylim(ax(7), [0 300]);
ylim(ax(8), [0 300]);

xlabel(ax(7:8), '$K_D$ (log$_{10}$ a.u.)', 'FontSize', 24)
ylabel(ax(7), 'EMax', 'FontSize', 24)
ylabel(ax(8), 'EMax', 'FontSize', 24)

% Axis Ticks and Format
set(ax(1:8),'GridLineWidth',0.5,'Linewidth',2)
xlbl = 10.^(-3:2:3);
xlabels = {'-3', '-1', '1', '3'};
xticks(ax(1:6),xlbl);
xticklabels(ax(1:6), xlabels)
xsecondarylabel([ax(2) ax(5)], '')
set(ax(1:6),'TickLabelInterpreter', 'Latex')
yticklabels(ax(2:3), '')
yticklabels(ax(5:6), '')

% Legend Format
l1 = legend(ax(7), [utop{:,1}], arrayfun(@num2str,kp,'UniformOutput',0),'Units','normalized',...
        'Position',[0.88 0.76 0.1 0.2], 'FontSize', 20);
l2 = legend(ax(8), [ubot{:,1}], arrayfun(@num2str,kp,'UniformOutput',0),'Units','normalized',...
        'Position',[0.88 0.52 0.1 0.2], 'FontSize', 20);
title(l1, '$k_p$');
title(l2, '$k_p$');

% Correct axes position
for i=1:6
    set(ax(i),'Units','normalized','Position',[left(i) bottom(i) width height]);
end

%set(findall(f,'-property','FontSize'),'FontSize',20)
set(findall(f,'-property','Interpreter'),'Interpreter','Latex')

figure(f);

% Save
set(f, 'Units', 'inches')
pos = get(f, 'Position');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(f, '../Figures/Manuscript/Final-Figures/KPR2','-dpdf','-r0');

%% Bound NP versus Bound TCR comparison

clear all
f = figure('Color', 'white','Position',[10 100 1500 1500]);

tlo = tiledlayout(1,2); 

ax(1) = nexttile(tlo); ax(2) = nexttile(tlo);

% New Simulations
koff_vals = [0.0001, 0.001, 0.01];
kd_vals = koff_vals ./ 0.1;
rho_vals = linspace(-4,3.5,16);
bound_type = {'phos', 'bound', 'np', 'covered','fracphos', 'normBound', 'normPhos','rescaled'};

tpc = 20;

valence = [1,2,3,4,5,10];

for kd = kd_vals
    TCR_NP_Comparison(ax(1:2), kd_vals, rho_vals, tpc, valence, colors)
end
%% Surface Capacity
clear all

f1 = figure('Color','white','Position',[10 100 1100 375]);
f2 = figure('Color','white','Position',[10 100 1100 375]);

left = [0.1, 0.52, 0.1, 0.52];
bottom = [0.2, 0.2, 0.2, 0.2];
width = 0.30;
height = 0.72;

for i=1:2
    ax(i) = axes(f1, 'Units','normalized','Position',[left(i) bottom(i) width height]);
end

for i=3:4
    ax(i) = axes(f2, 'Units','normalized','Position',[left(i) bottom(i) width height]);
end

tcrs_per_cluster = [1,3,5,10,20];
radius = [8,14,20];
col = abyss(length(radius));

hold(ax(:),'on');
grid(ax(:),'on');
for i = 1:length(radius)
    plotCapacity(ax(1), tcrs_per_cluster,radius(i), 'np', col(i,:), true)
    plotCapacity(ax(2), tcrs_per_cluster,radius(i), 'covered', col(i,:), true)
    plotCapacity(ax(3), tcrs_per_cluster,radius(i), 'bound', col(i,:), true)
    plotCapacity(ax(4), tcrs_per_cluster,radius(i), 'nc', col(i,:), true)
end

set(ax(1:4),'GridLineWidth',0.5,'Linewidth',2)
xticks(ax(:), tcrs_per_cluster);
l1 = legend(ax(2), 'Position',[0.84 0.66 0.15 0.15]);
l2 = legend(ax(3), 'Location','northeast');

title([l1 l2], 'NP Radius');

set(findall(f1,'-property','FontSize'),'FontSize',20)
set(findall(f1,'-property','Interpreter'),'Interpreter','Latex')
set(findall(f2,'-property','FontSize'),'FontSize',20)
set(findall(f2,'-property','Interpreter'),'Interpreter','Latex')

% Save
set(f1, 'Units', 'inches')
pos = get(f1, 'Position');
set(f1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
%print(f1, '../Figures/Manuscript/Final-Figures/CarryingCapacity','-dpdf','-r0');

%% Hill Coefficient Figure

f2 = figure('Color','white');
ax(1) = axes(f2);

%Plotting

col = gray(length(kd_vals));
hol = winter(length(kd_vals));
sol = autumn(length(kd_vals));

file = 'TauLeap';

hold(ax(1),'on')
HillCoeffPlots(ax(1), kd_vals, rho_vals, 5, 3, bound_type{6}, hol, file)
HillCoeffPlots(ax(1), kd_vals, rho_vals, 5, 1, bound_type{6}, col, file)
HillCoeffPlots(ax(1), kd_vals, rho_vals, 5, 20, bound_type{6}, sol, file)

set(ax(1), 'xscale','log','tickdir','out')
xlim(ax(1), [0.0005 2]);
ylim(ax(1), [0 2]);

set(findall(f2,'-property','FontSize'),'FontSize',16)
set(findall(f2,'-property','Interpreter'),'Interpreter','Latex')

%% Bound TCRs per NP

f3 = figure('Color', 'white');

ax(1) = subplot(2,4,1:2); ax(3) = subplot(2,4,3:4);
ax(2) = subplot(2,4,5:6); ax(4) = subplot(2,4,7:8);

hols = autumn(length(koff_vals)+1);
hols = hols(1:length(koff_vals),:);
cols = winter(length(koff_vals));

for i=1:4 
    % Axis scale for DR Curves
    set(ax(i),'xscale','log','tickdir','out')
    ylim(ax(i), [0 5])
end

scat1 = DoseResponsePlots(ax(1), kd_vals(1), [1,1], rho_vals, 1, 5, bound_type{4}, hols(1,:),true,true,20);
scat2 = DoseResponsePlots(ax(2), kd_vals(2), [1,1], rho_vals, 3, 5, bound_type{4}, hols(2,:), true,true,20);
scat3 = DoseResponsePlots(ax(3), kd_vals(3), [1,1], rho_vals, 20,5, bound_type{4}, hols(3,:),true,true,14);
scat4 = DoseResponsePlots(ax(4), kd_vals(4), [1,1], rho_vals,20,5, bound_type{4}, hols(4,:),true, true,8);

set(findall(f3,'-property','FontSize'),'FontSize',16)
set(findall(f3,'-property','Interpreter'),'Interpreter','Latex')

%% Hill Coefficient Plots

function y = HillCoeffPlots(ax, kd_vals, rho_vals, valence, tpc, bound_type, colors, fname)
    x_vals = linspace(min(rho_vals), max(rho_vals), 1000);
    rho_vals = 10.^rho_vals;
    x_vals = 10.^x_vals;
    
    if tpc == 1
        surface_type = 'uniform';
        file = [fname,'/0TPC/v',num2str(valence)];
    else
        surface_type = 'clusters';
        file = [fname,'/',num2str(tpc),'TPC/v',num2str(valence)];
    end
    
    hill_cof = [];
    
    if strcmp(bound_type,'bound')
        k=1;
    elseif strcmp(bound_type,'phos')
        k=2;
    else
        k=3;
    end
    
    i=1;
    
    for koff = kd_vals ./ 10
        [mean_data, var_data, output_fit] = load_DR(koff, rho_vals, surface_type, file, true);
        fit_data = output_fit{k};
        hill_cof = [hill_cof, fit_data.n];
    end
    
    x = log10(kd_vals); y = hill_cof';
    X = [ones(length(x),1),x'];

    b = X \ y;
    
    hold(ax,'on');
    scatter(ax, kd_vals, hill_cof, [], colors, 'filled');
    %plot(ax, kd_vals, log10(10.^b(1)*kd_vals.^b(2)));
    title(ax, ['Slope = ',num2str(b(2))]);
    ylabel(ax,'Hill Coefficient');
    xlabel(ax,'$K_D$');
    grid(ax,'on');
end


%% Fitting Function

function [meanC, meanU, varC, varU, c_ft, u_ft] = Hill_fitting(koff, rho_vals, surface_type)

    [meanC, meanU, varC, varU] = DR_fitting(koff, rho_vals, surface_type);
    

    % Fit Hill Function to dose response

    fo = fitoptions('Method','NonlinearLeastSquares',...
                   'Lower',[0,0,0],...
                   'Upper',[Inf,max(rho_vals),10],...
                   'StartPoint',[1 1 1]);
    ft = fittype('eMax*x^n / (x^n + ec50^n)','options',fo);

    [c_ft,c_gof] = fit(rho_vals',meanC(2,:)',ft);

    [u_ft,u_gof] = fit(rho_vals',meanU(2,:)',ft);

end