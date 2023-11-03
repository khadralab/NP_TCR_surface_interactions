clear all

tstart = tic;

% NP and Kinetic Params
rNP = 20;
vh = 5;
k0 = 0.1;
kon = 0.1;
koff = 0.5;
T0 = 1;
rho = 1;
np_params = [rNP, vh, k0, kon, koff, T0, rho];

% Simulation Params
total_t = 5000;
fps = 10;
sim_params = [total_t, fps];
total_sims = 500;
final_time = linspace(0,total_t,total_t*fps);

% Strong ligand: 1-10s --> koff = 0.1-1
% Weak ligand: 0.1-3s --> koff = 0.3-10

i=1;
cT = [];
uT = [];
rho_vals = linspace(0.5,4,8);
rho_vals = 10.^rho_vals;
Time = zeros(1,length(rho_vals));

for rho = rho_vals
    tic
    
     % Initialize variables
    cluster_bound_tcr = [];
    cluster_bound_np = [];
    cluster_phos_tcr = [];
    cluster_time = [];
    homo_bound_tcr = [];
    homo_bound_np = [];
    homo_phos_tcr = [];
    homo_time = [];
    np_params = [rNP, vh, k0, kon, koff, T0, rho];
    
    disp(['rho = ',num2str(rho)]);
    
    % Cluster Surface
    rSurf = 1000;
    num_clusters = 15;
    cluster_radius = 50;
    tcr_per_cluster = 20;
    num_tcr = 300;
    rTCR = 5;
    tcr_params = [rSurf, num_clusters, cluster_radius, tcr_per_cluster, num_tcr, rTCR];
    
    cfile = ['LongSims/Dose-Response/koff50/Cluster_rho',num2str(round(rho*100))];
    ufile = ['LongSims/Dose-Response/koff50/Uniform_rho',num2str(round(rho*100))];
    
    nsims=0;
    while nsims < total_sims
        [bound_tcr, bound_np, phos_tcr, states, time] = Gillespie_KPR(tcr_params, np_params, sim_params);
        bound_tcr = interp1(time, bound_tcr, final_time, 'previous');
        bound_np = interp1(time, bound_np, final_time, 'previous');
        phos_tcr = interp1(time, phos_tcr, final_time, 'previous');
        cluster_bound_tcr = [cluster_bound_tcr; bound_tcr];
        cluster_bound_np = [cluster_bound_np; bound_np];
        cluster_phos_tcr = [cluster_phos_tcr; phos_tcr];
        nsims = nsims+1;
    end
        
    % Homogeneous Surface
    rSurf = 1000;
    num_clusters = 0;
    cluster_radius = 50;
    tcr_per_cluster = 20;
    num_tcr = 300;
    rTCR = 5;
    tcr_params = [rSurf, num_clusters, cluster_radius, tcr_per_cluster, num_tcr, rTCR];
    file = ['Homo_v',num2str(vh),'r',num2str(rNP),'koff',num2str(koff*100)];

    nsims = 0;
    while nsims < total_sims
        [bound_tcr, bound_np, phos_tcr, homo_states, time] = Gillespie_KPR(tcr_params, np_params, sim_params);
        bound_tcr = interp1(time, bound_tcr, final_time, 'previous');
        bound_np = interp1(time, bound_np, final_time, 'previous');
        phos_tcr = interp1(time, phos_tcr, final_time, 'previous');
        homo_bound_tcr = [homo_bound_tcr; bound_tcr];
        homo_bound_np = [homo_bound_np; bound_np];
        homo_phos_tcr = [homo_phos_tcr; phos_tcr];
        nsims = nsims + 1;
    end
    
    save(ufile,'homo_bound_tcr','homo_bound_np','homo_phos_tcr')
    save(cfile,'cluster_bound_tcr','cluster_bound_np','cluster_phos_tcr')
    
    Time(i) = toc
    i=i+1;
    
end

tend = toc(tstart)
