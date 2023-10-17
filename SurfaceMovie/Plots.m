%% Clusters versus Uniform (side-by-side)

v=2; r=20; koff=0.05; rho=1;

load(['LongSims/Cluster_v',num2str(v),'rho',num2str(rho*10),'koff',num2str(koff*100),'.mat'])
load(['LongSims/Homo_v',num2str(v),'rho',num2str(rho*10),'koff',num2str(koff*100),'.mat'])

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
plot(final_time, cluster_bound_np,'-','Color', [0, 0.5, 0.8, 0.01])
hold on
plot(final_time, mean(cluster_bound_np,1), '-','Color', [0, 0.5, 0.8, 1],'Linewidth',2,'DisplayName','Clusters')
plot(final_time, homo_bound_np, '-','Color', [0.8, 0, 0.5, 0.01])
plot(final_time, mean(homo_bound_np,1), '-','Color', [0.8, 0, 0.5, 1],'Linewidth',2,'DisplayName','Uniform')
xlim([0 5000])
xlabel('Time')
ylabel('Bound NPs')

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
histogram(homo_bound_np, edges, 'Normalization','pdf','FaceColor', [0.8, 0, 0.5], 'FaceAlpha', 0.5,'DisplayName','Uniform')
hold on
histogram(cluster_bound_np, edges, 'Normalization','pdf','FaceColor', [0, 0.5, 0.8],'FaceAlpha', 0.5,'DisplayName','Clusters')
ylabel('PDF')
xlabel('Bound NPs')
xlim([0 10])

legend()
set(findall(f1,'-property','FontSize'),'FontSize',16)


f2 = figure(2);
subplot(221)
plot(final_time, cluster_phos_tcr,'-','Color', [0, 0.5, 0.8, 0.01])
hold on
plot(final_time, homo_phos_tcr,'-','Color', [0.8, 0, 0.5, 0.01])
plot(final_time, mean(cluster_phos_tcr,1), '-','Color', [0, 0.5, 0.8, 1],'Linewidth',2,'DisplayName','Clusters')
plot(final_time, mean(homo_phos_tcr,1), '-','Color', [0.8, 0, 0.5, 1],'Linewidth',2,'DisplayName','Uniform')
xlim([0 5000])
xlabel('Time')
ylabel('Phos. TCRs')

subplot(222)
histogram(homo_phos_tcr, edges, 'Normalization','pdf','FaceColor', [0.8, 0, 0.5], 'FaceAlpha', 0.5,'DisplayName','Uniform')
hold on
histogram(cluster_phos_tcr, edges, 'Normalization','pdf','FaceColor', [0, 0.5, 0.8],'FaceAlpha', 0.5,'DisplayName','Clusters')
ylabel('PDF')
xlabel('Phos. TCRs')
xlim([0 20])

legend()
set(findall(f2,'-property','FontSize'),'FontSize',16)

%% Specificity Plots (strong versus weak ligands)

cols = [[0.3, 0.8, 0.5]; [0.8, 0.4, 0]];
i=1;

for koff = [0.05, 0.5]
    load(['LongSims/Cluster_v',num2str(v),'rho',num2str(rho*10),'koff',num2str(koff*100),'.mat'])
    load(['LongSims/Homo_v',num2str(v),'rho',num2str(rho*10),'koff',num2str(koff*100),'.mat'])
    
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
    histogram(cluster_bound_np, edges, 'Normalization','pdf','FaceColor', cols(i,:))
    hold on
    xlim([0 10])
    ylim([0 1])
    xlabel('Bound NP')
    
    subplot(224)
    histogram(homo_bound_np, edges, 'Normalization','pdf','FaceColor', cols(i,:))
    hold on
    xlim([0 10])
    ylim([0 1])
    xlabel('Bound NP')
    
    legend('Strong', 'Weak')
    
    f4 = figure(4);
    subplot(221)
    histogram(cluster_phos_tcr, edges, 'Normalization','pdf','FaceColor', cols(i,:))
    hold on
    ylabel('PDF')
    xlabel('Phos. TCRs')
    title('Clusters')
    xlim([0 20])

    subplot(222)
    histogram(homo_phos_tcr, edges, 'Normalization','pdf','FaceColor', cols(i,:))
    hold on
    ylabel('PDF')
    xlabel('Phos. TCRs')
    title('Uniform')
    xlim([0 20])
    
    i=i+1;
end

legend('Strong', 'Weak')
set(findall(f3,'-property','FontSize'),'FontSize',16)
set(findall(f4,'-property','FontSize'),'FontSize',16)

%%

cols = [[0.8, 0.5, 0.5]; [0., 0.5, 0.8]];

mean_bound = zeros(3,2);
dev_bound = zeros(3,2);

i=1;

rho = 1;
koff = 0.5;

for v = [1,2,5]
    load(['LongSims/Cluster_v',num2str(v),'rho',num2str(rho*10),'koff',num2str(koff*100),'.mat'])
    load(['LongSims/Homo_v',num2str(v),'rho',num2str(rho*10),'koff',num2str(koff*100),'.mat'])
    
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

