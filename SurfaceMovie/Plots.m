%% Clusters versus Uniform (side-by-side)

v=2; r=20; koff=0.05; rho=1;

load(['LongSims/Original/Cluster_v',num2str(v),'rho',num2str(rho*10),'koff',num2str(koff*100),'.mat'])
load(['LongSims/Original/Homo_v',num2str(v),'rho',num2str(rho*10),'koff',num2str(koff*100),'.mat'])

% Simulation Params
total_t = 5000;
fps = 10;
final_time = linspace(0,total_t,total_t*fps);

% Time Series
f1 = figure(1);

subplot(221)
plot(final_time, cluster_bound_tcr,'-','Color', [0, 0.5, 0.8, 0.01])
hold on
plot(final_time, homo_bound_tcr,'-','Color', [0.8, 0, 0.5, 0.01])
plot(final_time, mean(cluster_bound_tcr,1), '-','Color', [0, 0.5, 0.8, 1],'Linewidth',2,'DisplayName','Clusters')
plot(final_time, mean(homo_bound_tcr,1), '-','Color', [0.8, 0, 0.5, 1],'Linewidth',2,'DisplayName','Uniform')
xlim([0 5000])
xlabel('Time')
ylabel('Bound TCRs')

subplot(222)
plot(final_time, cluster_phos_tcr,'-','Color', [0, 0.5, 0.8, 0.01])
hold on
plot(final_time, homo_phos_tcr,'-','Color', [0.8, 0, 0.5, 0.01])
plot(final_time, mean(cluster_phos_tcr,1), '-','Color', [0, 0.5, 0.8, 1],'Linewidth',2,'DisplayName','Clusters')
plot(final_time, mean(homo_phos_tcr,1), '-','Color', [0.8, 0, 0.5, 1],'Linewidth',2,'DisplayName','Uniform')
xlim([0 5000])
xlabel('Time')
ylabel('Phos. TCRs')

% Histograms

edges = linspace(0, 20, 21) - 0.5;

subplot(223)
histogram(homo_bound_tcr, edges, 'Normalization','pdf','FaceColor', [0.8, 0, 0.5], 'FaceAlpha', 0.5)
hold on
histogram(cluster_bound_tcr, edges, 'Normalization','pdf','FaceColor', [0, 0.5, 0.8],'FaceAlpha', 0.5)
ylabel('PDF')
xlabel('Bound TCR')
xlim([0 10])

subplot(224)
histogram(homo_phos_tcr, edges, 'Normalization','pdf','FaceColor', [0.8, 0, 0.5], 'FaceAlpha', 0.5,'DisplayName','Uniform')
hold on
histogram(cluster_phos_tcr, edges, 'Normalization','pdf','FaceColor', [0, 0.5, 0.8],'FaceAlpha', 0.5,'DisplayName','Clusters')
ylabel('PDF')
xlabel('Phos. TCRs')
xlim([0 20])

legend()
set(findall(f1,'-property','FontSize'),'FontSize',16)

%% Specificity Plots (strong versus weak ligands)

cols = [[0.3, 0.8, 0.5]; [0.8, 0.4, 0]];
i=1;

for koff = [0.05, 0.5]
    load(['LongSims/Original/Cluster_v',num2str(v),'rho',num2str(rho*10),'koff',num2str(koff*100),'.mat'])
    load(['LongSims/Original/Homo_v',num2str(v),'rho',num2str(rho*10),'koff',num2str(koff*100),'.mat'])
    
    f3 = figure(3);
    subplot(221)
    histogram(cluster_bound_tcr, edges, 'Normalization','pdf','FaceColor', cols(i,:))
    hold on
    xlabel('Bound TCRs')
    title('Clusters')
    xlim([0 10])
    ylim([0 1])
    
    subplot(222)
    histogram(homo_bound_tcr, edges, 'Normalization','pdf','FaceColor', cols(i,:))
    hold on
    xlim([0 10])
    ylim([0 1])
    xlabel('Bound TCRs')
    title('Uniform')
    
    subplot(223)
    histogram(cluster_phos_tcr, edges, 'Normalization','pdf','FaceColor', cols(i,:))
    hold on
    ylabel('PDF')
    xlabel('Phos. TCRs')
    title('Clusters')
    xlim([0 10])

    subplot(224)
    histogram(homo_phos_tcr, edges, 'Normalization','pdf','FaceColor', cols(i,:))
    hold on
    ylabel('PDF')
    xlabel('Phos. TCRs')
    title('Uniform')
    xlim([0 10])
    
    legend('Strong', 'Weak')
    
    i=i+1;
end

legend('Strong', 'Weak')
set(findall(f3,'-property','FontSize'),'FontSize',16)
set(findall(f4,'-property','FontSize'),'FontSize',16)

%% Sensitivity: Valence

cols = [[0.8, 0.5, 0.5]; [0., 0.5, 0.8]];

mean_bound = zeros(3,2);
dev_bound = zeros(3,2);

i=1;

rho = 1;
koff = 0.5;

for v = [1,2,5]
    load(['LongSims/Original/Cluster_v',num2str(v),'rho',num2str(rho*10),'koff',num2str(koff*100),'.mat'])
    load(['LongSims/Original/Homo_v',num2str(v),'rho',num2str(rho*10),'koff',num2str(koff*100),'.mat'])
    
    mean_bound(i,1) = mean(cluster_bound_tcr(:,end),1);
    mean_bound(i,2) = mean(homo_bound_tcr(:,end),1);
    dev_bound(i,1) = std(cluster_bound_tcr(:,end),1);
    dev_bound(i,2) = std(homo_bound_tcr(:,end),1);
    
    i=i+1;
end


f5 = figure(5);
hb = bar(mean_bound);
xlabel('NP Valence')
ylabel('Mean Bound TCRs')
xticks([1 2 3])
xticklabels({'1','2','5'})
%yticks([0:2:40])
%yticklabels(num2cell([0:2:40]))

grid on
grid minor

hold on

for k = 1:2
    
    xpos = hb(k).XData + hb(k).XOffset;

    er = errorbar(xpos,mean_bound(:,k),dev_bound(:,k)/sqrt(500), 'LineStyle','none','Color',zeros(1,3),'LineWidth',1);    
    %er.LineStyle = 'none';
end

hold off
legend('Clusters','Uniform','Location','Northwest')

%ylim([0 38])
%breakyaxis([3,36])
set(findall(f5,'-property','FontSize'),'FontSize',16)

%% Sensitivity: Phosphorylation rate

cols = [[0.8, 0.5, 0.5]; [0., 0.5, 0.8]];

mean_bound = zeros(3,2);
dev_bound = zeros(3,2);

i=1;

rho = 1;
r=20;
v=5;

for kp = [1,5,10]
    load(['LongSims/SensitivityAnalysis/kp',num2str(kp),'/Cluster_v',num2str(v),'rho',num2str(rho*10),'kp',num2str(kp),'.mat'])
    load(['LongSims/SensitivityAnalysis/kp',num2str(kp),'/Homo_v',num2str(v),'rho',num2str(rho*10),'kp',num2str(kp),'.mat'])
    
    mean_bound(i,1) = mean(cluster_phos_tcr(:,end),1);
    mean_bound(i,2) = mean(homo_phos_tcr(:,end),1);
    dev_bound(i,1) = std(cluster_phos_tcr(:,end),1);
    dev_bound(i,2) = std(homo_phos_tcr(:,end),1);
    
    i=i+1;
end


f6 = figure(6);
hb = bar(mean_bound);
xlabel('Phosphorylation Rate')
ylabel('Mean Phosphorylated TCRs')
xticks([1 2 3])
xticklabels({'0.01','0.05','0.1'})
%yticks([0:2:40])
%yticklabels(num2cell([0:2:40]))

grid on
grid minor

hold on

for k = 1:2
    
    xpos = hb(k).XData + hb(k).XOffset;

    er = errorbar(xpos,mean_bound(:,k),dev_bound(:,k)/sqrt(500), 'LineStyle','none','Color',zeros(1,3),'LineWidth',1);    
    %er.LineStyle = 'none';
end

hold off
legend('Clusters','Uniform','Location','Northwest')

%ylim([0 38])
%breakyaxis([3,36])
set(findall(f6,'-property','FontSize'),'FontSize',16)

%% Sensitivity: TCR Cross-linking rate

cols = [[0.8, 0.5, 0.5]; [0., 0.5, 0.8]];

rho_vals = linspace(-2,0.5,6);
rho_vals = 10.^rho_vals;

% k0 = (0.01, 0.05, 0.1, 0.5) * 100

k0_vals = [1, 5];

mean_bound = zeros(length(k0_vals)*2,length(rho_vals));
dev_bound = zeros(length(k0_vals)*2,length(rho_vals));

i=1;

v=5;

for k0 = k0_vals
    for j = [1:length(rho_vals)]
        rho = rho_vals(j);
        
        load(['LongSims/SensitivityAnalysis/k0',num2str(k0),'/Cluster_v',num2str(v),'rho',num2str(round(rho*100)),'k0',num2str(k0),'.mat'])
        load(['LongSims/SensitivityAnalysis/k0',num2str(k0),'/Homo_v',num2str(v),'rho',num2str(round(rho*100)),'k0',num2str(k0),'.mat'])

        mean_bound(i,j) = mean(cluster_bound_tcr(:,end),1);
        mean_bound(i+1,j) = mean(homo_bound_tcr(:,end),1);
        %dev_bound(i,1) = std(cluster_bound_tcr(:,end),1);
        %dev_bound(i,2) = std(homo_bound_tcr(:,end),1);
    end
    
    i=i+2;
end


f7 = figure(7);
semilogx(rho_vals, mean_bound)
xlabel('NP Dose')
ylabel('Mean Bound TCRs')
%yticks([0:2:40])
%yticklabels(num2cell([0:2:40]))
grid on
grid minor
hold off
legend('C1', 'U1', 'C5', 'U5')

%ylim([0 38])
%breakyaxis([3,36])
set(findall(f7,'-property','FontSize'),'FontSize',16)

%% Dose-Response Curves
% Plots of results from DoseResponse.m for: (r,v, koff) = (20,5,0.05)

rho_vals = linspace(-2,1,7);
rho_vals = 10.^rho_vals;

hols = autumn(3);
cols = winter(3);

f8 = figure(8);
i=1;

for koff = [0.05, 0.5]
    [aC, aU, c3, u3] = DR_fitting(koff, rho_vals);

    subplot(211)
    grid on
    semilogx(rho_vals, aC, 'Color', hols(i,:), 'DisplayName',['Cluster ',num2str(koff)],'Linewidth',1)
    hold on
    semilogx(rho_vals, aU, 'Color', cols(i,:), 'DisplayName',['Uniform ',num2str(koff)],'Linewidth',1)
    
    subplot(212)
    grid on
    loglog(rho_vals, aC, 'Color', hols(i,:), 'DisplayName',['Cluster ',num2str(koff)],'Linewidth',1)
    hold on
    loglog(rho_vals, aU, 'Color', cols(i,:), 'DisplayName',['Uniform ',num2str(koff)],'Linewidth',1)
    
    i=i+1;
end

hold off
xlabel('NP Concentration')
ylabel('Bound TCRs')
legend()

set(findall(f8,'-property','FontSize'),'FontSize',16)

%% Heatmap Dose Response????

rho_vals = [0.01, 0.03, 0.1, 0.31, 1];
koff_vals = [0.05, 0.5];

cluster_map = zeros(length(koff_vals), length(rho_vals));
uniform_map = zeros(length(koff_vals), length(rho_vals));

for i = 1:length(koff_vals)
    koff = koff_vals(i);
    
    [aC, aU, c3, u3] = DR_fitting(koff, rho_vals);
    
    cluster_map(i,:) = aC;
    uniform_map(i,:) = aU;
end

f9 = figure(9);

subplot(211)
imagesc(koff_vals, log10(rho_vals), cluster_map)

subplot(212)
imagesc(koff_vals, log10(rho_vals), uniform_map)
xlabel('NP Concentration')
ylabel('Dissociation Constant (Koff)')


%% Local Functions

function [aC, aU, c3, u3] = DR_fitting(koff, rho_vals)
    
    aC = []; vC = [];
    aU = []; vU = [];

    for rho = rho_vals
        load(['LongSims/Dose-Response/koff',num2str(koff*100),'/Cluster_rho',num2str(floor(rho*100)),'.mat'])

        aC = [aC, mean(cluster_bound_tcr(:,end))];
        vC = [vC, std(cluster_bound_tcr(:,end))];

        load(['LongSims/Dose-Response/koff',num2str(koff*100),'/Uniform_rho',num2str(floor(rho*100)),'.mat'])

        aU = [aU, mean(homo_bound_tcr(:,end))];
        vU = [vU, std(homo_bound_tcr(:,end))];

    end

    % Fit Hill Function to dose response

    fo = fitoptions('Method','NonlinearLeastSquares',...
                   'Lower',[0,0],...
                   'Upper',[Inf,max(rho_vals)],...
                   'StartPoint',[1 1]);
    ft = fittype('eMax*x^n / (x^n + ec50^n)','problem','n','options',fo);

    [c2,gof2] = fit(rho_vals',aC',ft,'problem',2);
    [c3,gof2] = fit(rho_vals',aC',ft,'problem',1);

    [u2,gof2] = fit(rho_vals',aU',ft,'problem',2);
    [u3,gof2] = fit(rho_vals',aU',ft,'problem',1);

    
end


