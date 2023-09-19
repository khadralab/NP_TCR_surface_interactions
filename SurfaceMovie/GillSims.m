%% Multiple Runs

clear all

% Surface Params
rSurf = 1000;
num_clusters = 15;
cluster_radius = 50;
tcr_per_cluster = 20;
num_tcr = 300;
rTCR = 5;

% NP and Kinetic Params
rNP = 20;
vh = 5;
k0 = 0.1;
kon = 0.1;
koff = 0.05;
T0 = 1;
rho = 1;
np_params = [rNP, vh, k0, kon, koff, T0, rho];

% Simulation Params
total_t = 5000;
fps = 10;
sim_params = [total_t, fps];
total_sims = 100;
final_time = linspace(0,total_t,total_t*fps);

% Strong ligand: 1-10s --> koff = 0.1-1
% Weak ligand: 0.1-3s --> koff = 0.3-10

i=1;
for num_clusters = [15,0]
    
    % Initialize variables
    cluster_bound_tcr = [];
    cluster_bound_np = [];
    cluster_phos_tcr = [];
    cluster_time = [];
    homo_bound_tcr = [];
    homo_bound_np = [];
    homo_phos_tcr = [];
    homo_time = [];
        
    % Strong Ligand
    koff = 0.05;
    np_params = [rNP, vh, k0, kon, koff, T0, rho];
    tcr_params = [rSurf, num_clusters, cluster_radius, tcr_per_cluster, num_tcr, rTCR];    
    nsims=0;
    
    while nsims < total_sims
        [bound_tcr, bound_np, phos_tcr, time] = Gillespie_KPR(tcr_params, np_params, sim_params);
        bound_tcr = interp1(time, bound_tcr, final_time, 'previous');
        bound_np = interp1(time, bound_np, final_time, 'previous');
        phos_tcr = interp1(time, phos_tcr, final_time, 'previous');
        cluster_bound_tcr = [cluster_bound_tcr; bound_tcr];
        cluster_bound_np = [cluster_bound_np; bound_np];
        cluster_phos_tcr = [cluster_phos_tcr; phos_tcr];
        nsims = nsims+1;
    end
    
    % Weak Ligand
    koff = 0.5;
    np_params = [rNP, vh, k0, kon, koff, T0, rho];
    nsims = 0;
    
    while nsims < total_sims
        [bound_tcr, bound_np, phos_tcr, time] = Gillespie_KPR(tcr_params, np_params, sim_params);
        bound_tcr = interp1(time, bound_tcr, final_time, 'previous');
        bound_np = interp1(time, bound_np, final_time, 'previous');
        phos_tcr = interp1(time, phos_tcr, final_time, 'previous');
        homo_bound_tcr = [homo_bound_tcr; bound_tcr];
        homo_bound_np = [homo_bound_np; bound_np];
        homo_phos_tcr = [homo_phos_tcr; phos_tcr];
        nsims = nsims + 1;
    end
    
    figure()
    if i == 1
        sgtitle(['Clusters: r=',num2str(rNP),', v=',num2str(vh)])
        file = ['Clusters_koff'];
    else
        sgtitle(['Uniform: r=',num2str(rNP),', v=',num2str(vh)])
        file = ['Uniform_koff'];
    end
    
    subplot(221)
    histogram(homo_bound_tcr,'Normalization','pdf','FaceColor', [0.8, 0.4, 0], 'FaceAlpha', 0.5)
    hold on
    histogram(cluster_bound_tcr,'Normalization','pdf','FaceColor', [0.3, 0.8, 0.5],'FaceAlpha', 0.5)
    ylabel('pdf')
    title('Bound TCR')

    subplot(222)
    histogram(homo_bound_np,'Normalization','pdf','FaceColor', [0.8, 0.4, 0], 'FaceAlpha', 0.5)
    hold on
    histogram(cluster_bound_np,'Normalization','pdf','FaceColor', [0.3, 0.8, 0.5],'FaceAlpha', 0.5)
    ylabel('pdf')
    title('Bound NPs')

    subplot(223)
    histogram(homo_phos_tcr,'Normalization','pdf','FaceColor', [0.8, 0.4, 0], 'FaceAlpha', 0.5, 'DisplayName', 'Weak Ligand')
    hold on
    histogram(cluster_phos_tcr,'Normalization','pdf','FaceColor', [0.3, 0.8, 0.5], 'FaceAlpha', 0.5, 'DisplayName','Strong Ligand')
    ylabel('pdf')
    title('Phosph. TCRs')
    legend(gca,'Position',[0.5 0.35 0.2 0.1])
    
    %savefig(gcf,file);
    %saveas(gcf,file,'png');
    
    i=i+1;
end