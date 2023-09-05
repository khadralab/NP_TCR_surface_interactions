clear all

rNP = 50;
vh = 5;
k0 = 0.01;
kon = 0.01;
koff = 0.05;
T0 = 1;
rho = 1;

np_params = [rNP, vh, k0, kon, koff, T0, rho];

% Clustered Surface
rSurf = 5000;
num_clusters = 15;
cluster_radius = 50;
tcr_per_cluster = 20;
num_tcr = 300;
rTCR = 5;

tcr_params = [rSurf, num_clusters, cluster_radius, tcr_per_cluster, num_tcr, rTCR];

% Sim and Save Params
total_t = 5000;
fps = 10;
file = 'mov4';

sim_params = [total_t, fps];

[bound_tcr, bound_np, phos_tcr, time] = Gillespie_KPR(tcr_params, np_params, sim_params, file);

% Homogeneous Surface
rSurf = 1000;
num_clusters = 0;
cluster_radius = 50;
tcr_per_cluster = 20;
num_tcr = 300;
rTCR = 5;

tcr_params = [rSurf, num_clusters, cluster_radius, tcr_per_cluster, num_tcr, rTCR];

[bound_tcr2, bound_np2, phos_tcr2, time2] = Gillespie_KPR(tcr_params, np_params, sim_params, file);

figure()
subplot(221)
plot(time, bound_np)
hold on
plot(time2, bound_np2)
xlabel('Time')
title('Bound NPs')

subplot(222)
plot(time, bound_tcr)
hold on
plot(time2, bound_tcr2)
xlabel('Time')
title('Bound TCRs')

subplot(223)
plot(time,phos_tcr)
hold on
plot(time2, phos_tcr2)
xlabel('Time')
title('Phosphorylated TCRs')

subplot(224)
plot(time, cumsum(phos_tcr),'DisplayName','Clusters')
hold on
plot(time2, cumsum(phos_tcr2),'DisplayName','Homogeneous')
xlabel('Time')
title('Cumulative Phos. TCRs')
legend()


%% Multiple Runs

clear all

% NP and Kinetic Params
rNP = 20;
vh = 5;
k0 = 0.1;
kon = 0.1;
koff = 5;
T0 = 1;
rho = 1;

np_params = [rNP, vh, k0, kon, koff, T0, rho];

% Simulation Params

total_t = 500;
fps = 10;

sim_params = [total_t, fps];
total_sims = 100;

final_time = linspace(0,total_t,total_t*fps);

% Strong ligand: 1-10s --> koff = 0.1-1
% Weak ligand: 0.1-3s --> koff = 0.3-10

i=0;
for koff = [0.05,0.1,0.5,1]
    
    % Initialize variables
    cluster_bound_tcr = [];
    cluster_bound_np = [];
    cluster_phos_tcr = [];
    cluster_time = [];
    homo_bound_tcr = [];
    homo_bound_np = [];
    homo_phos_tcr = [];
    homo_time = [];
    
    disp(['Koff = ',num2str(koff)]);
    
    np_params = [rNP, vh, k0, kon, koff, T0, rho];
    
    % Cluster Surface
    rSurf = 1000;
    num_clusters = 15;
    cluster_radius = 50;
    tcr_per_cluster = 20;
    num_tcr = 300;
    rTCR = 5;

    tcr_params = [rSurf, num_clusters, cluster_radius, tcr_per_cluster, num_tcr, rTCR];
    
    file = ['Cluster_v',num2str(vh),'r',num2str(rNP),'koff',num2str(koff*100)];
    
    nsims=0;
    
    while nsims < total_sims
        [bound_tcr, bound_np, phos_tcr, time] = Gillespie_KPR(tcr_params, np_params, sim_params, file);
        
        bound_tcr = interp1(time, bound_tcr, final_time, 'previous');
        bound_np = interp1(time, bound_np, final_time, 'previous');
        phos_tcr = interp1(time, phos_tcr, final_time, 'previous');
        
        cluster_bound_tcr = [cluster_bound_tcr; bound_tcr];
        cluster_bound_np = [cluster_bound_np; bound_np];
        cluster_phos_tcr = [cluster_phos_tcr; phos_tcr];
        
        nsims = nsims+1;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    f1 = figure(2*i+1);
    subplot(221)
    hold on
    plot(final_time, cluster_bound_tcr,'-','Color', [0, 0.5, 0.8, 0.1])
    plot(final_time, mean(cluster_bound_tcr),'-','Color', [0, 1, 1, 1],'Linewidth',2)
    xlabel('Time')
    title('Bound TCR')

    subplot(222)
    hold on
    plot(final_time, cluster_bound_np,'-','Color', [0, 0.5, 0.8, 0.1])
    plot(final_time, mean(cluster_bound_np),'-','Color', [0, 1, 1, 1],'Linewidth',2)
    xlabel('Time')
    title('Bound NPs')

    subplot(223)
    hold on
    plot(final_time, cluster_phos_tcr,'-','Color', [0, 0.5, 0.8, 0.1])
    plot(final_time, mean(cluster_phos_tcr),'-','Color', [0, 1, 1, 1],'Linewidth',2)
    xlabel('Time')
    title('Phosph. TCRs')
    
    subplot(224)
    hold on
    plot(final_time, cumsum(cluster_phos_tcr,2),'-','Color', [0, 0.5, 0.8, 0.1])
    plot(final_time, mean(cumsum(cluster_phos_tcr,2)),'-','Color', [0, 1, 1, 1],'Linewidth',2)
    xlabel('Time')
    title('Sum Phosph. TCRs')
    
    f2 = figure(2*i+2);
    subplot(221)
    histogram(cluster_bound_tcr,'Normalization','pdf','FaceColor', [0, 0.5, 0.8],'FaceAlpha', 0.1)
    hold on
    title('Bound TCR')

    subplot(222)
    histogram(cluster_bound_np,'Normalization','pdf','FaceColor', [0, 0.5, 0.8],'FaceAlpha', 0.1)
    hold on
    title('Bound NPs')

    subplot(223)
    histogram(cluster_phos_tcr,'Normalization','pdf','FaceColor', [0, 0.5, 0.8], 'FaceAlpha', 0.1)
    hold on
    title('Phosph. TCRs')

    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
     % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
     
    % Homogeneous Surface
    rSurf = 1000;
    num_clusters = 0;
    cluster_radius = 50;
    tcr_per_cluster = 20;
    num_tcr = 300;
    rTCR = 5;
    
    file = ['Homo_v',num2str(vh),'r',num2str(rNP),'koff',num2str(koff*100)];

    tcr_params = [rSurf, num_clusters, cluster_radius, tcr_per_cluster, num_tcr, rTCR];
    
    nsims = 0;
    
    while nsims < total_sims
        [bound_tcr, bound_np, phos_tcr, time] = Gillespie_KPR(tcr_params, np_params, sim_params, file);
        
        bound_tcr = interp1(time, bound_tcr, final_time, 'previous');
        bound_np = interp1(time, bound_np, final_time, 'previous');
        phos_tcr = interp1(time, phos_tcr, final_time, 'previous');
        
        homo_bound_tcr = [homo_bound_tcr; bound_tcr];
        homo_bound_np = [homo_bound_np; bound_np];
        homo_phos_tcr = [homo_phos_tcr; phos_tcr];
        
        nsims = nsims + 1;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(f1)
    
    subplot(221)
    hold on
    plot(final_time, homo_bound_tcr,'-','Color', [0.8, 0, 0.5, 0.1])
    plot(final_time, mean(homo_bound_tcr),'-','Color', [1, 0, 1, 1],'Linewidth',2)
    hold off
    xlabel('Time')
    ylabel('Bound TCR')

    subplot(222)
    hold on
    plot(final_time, homo_bound_np,'-','Color', [0.8, 0, 0.5, 0.1])
    plot(final_time, mean(homo_bound_np),'-','Color', [1, 0, 1, 1],'Linewidth',2)
    hold off
    xlabel('Time')
    ylabel('Bound NPs')
    xlim([0 total_t])

    subplot(223)
    hold on
    plot(final_time, homo_phos_tcr,'-','Color', [0.8, 0, 0.5, 0.1])
    plot(final_time, mean(homo_phos_tcr),'-','Color', [1, 0, 1, 1],'Linewidth',2)
    hold off
    xlabel('Time')
    title('Phosph. TCRs')
    
    subplot(224)
    hold on
    plot(final_time, cumsum(homo_phos_tcr,2),'-','Color', [0.8, 0, 0.5, 0.1])
    plot(final_time, mean(cumsum(homo_phos_tcr,2)),'-','Color', [1, 0, 1, 1],'Linewidth',2)
    hold off
    xlabel('Time')
    title('Sum Phosph. TCRs')
    
    figure(f2)
    
    subplot(221)
    histogram(homo_bound_tcr,'Normalization','pdf','FaceColor', [0.8, 0, 0.5], 'FaceAlpha', 0.1)
    ylabel('pdf')
    title('Bound TCR')

    subplot(222)
    histogram(homo_bound_np,'Normalization','pdf','FaceColor', [0.8, 0, 0.5], 'FaceAlpha', 0.1)
    ylabel('pdf')
    title('Bound NPs')

    subplot(223)
    histogram(homo_phos_tcr,'Normalization','pdf','FaceColor', [0.8, 0, 0.5], 'FaceAlpha', 0.1)
    ylabel('pdf')
    title('Phosph. TCRs')
    
    sgtitle(f1, ['Koff = ',num2str(koff),', v = ',num2str(vh),', r = ',num2str(rNP)]);
    
    file = ['Compare_v',num2str(vh),'r',num2str(rNP),'koff',num2str(koff*100)];
    savefig(f1,file);
    saveas(f1,file,'png');
    
    %sgtitle(f2, [file(1:4),': Koff = ',num2str(koff),', v = ',num2str(vh),', r = ',num2str(rNP)]);
    file = ['Hist_v',num2str(vh),'r',num2str(rNP),'koff',num2str(koff*100)];
    saveas(f2,file,'png');
    
    i=i+1;
    
end



