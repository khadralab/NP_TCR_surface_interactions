function [cluster_states] = DoseResponseFunction(rho_index)
    
    % NP and Kinetic Params
    % Strong ligand: 1-10s --> koff = 0.1-1
    % Weak ligand: 0.1-3s --> koff = 0.3-10
    rNP = 20;
    vh = 5;
    k0 = 0.1;
    kon = 0.1;
    koff = 0.01;
    T0 = 1;

    % Simulation Params
    if rho_index < 0
        total_t = 1e5;
        tau = 1;
        sim_params = [total_t, tau];
    else
        total_t = 1e2;
        tau = 500;
        sim_params = [total_t, tau];
    end
    total_sims = 10;
    
    % NP concentrations
    rho_vals = linspace(-2,5,15);
    rho_vals = 10.^rho_vals;
    rho = rho_vals(rho_index);
    
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

    cfile = ['LongSims/temp/koff',num2str(koff*100),'/Cluster_rho',num2str(floor(rho*100))];
    ufile = ['LongSims/temp/koff',num2str(koff*100),'/Uniform_rho',num2str(floor(rho*100))];
    
    for nsims = [1:total_sims]
        [bound_tcr, bound_np, phos_tcr, cluster_states, time] = TauLeaping(tcr_params, np_params, sim_params);
        [time,j,~] = unique(time);
        total_t = floor(time(end))+1;
        final_time = linspace(0,total_t,total_t*10);

        bound_tcr = bound_tcr(j); bound_np = bound_np(j); phos_tcr = phos_tcr(j);
        
        bound_tcr = interp1(time, bound_tcr, final_time, 'previous');
        bound_np = interp1(time, bound_np, final_time, 'previous');
        phos_tcr = interp1(time, phos_tcr, final_time, 'previous');
        
        cluster_bound_tcr{nsims} = [bound_tcr; final_time];
        cluster_bound_np{nsims} = [bound_np; final_time];
        cluster_phos_tcr{nsims} = [phos_tcr; final_time];
    end
    
    save(cfile,'cluster_bound_tcr','cluster_bound_np','cluster_phos_tcr');
%{
    % Homogeneous Surface
    rSurf = 1000;
    num_clusters = 0;
    cluster_radius = 50;
    tcr_per_cluster = 20;
    num_tcr = 300;
    rTCR = 5;
    tcr_params = [rSurf, num_clusters, cluster_radius, tcr_per_cluster, num_tcr, rTCR];
    
    parfor nsims = [1:total_sims]
        [bound_tcr, bound_np, phos_tcr, uniform_states, time] = TauLeaping(tcr_params, np_params, sim_params);
        [time,j,~] = unique(time);
        total_t = floor(time(end))+1;
        final_time = linspace(0,total_t,total_t/tau);

        bound_tcr = bound_tcr(j); bound_np = bound_np(j); phos_tcr = phos_tcr(j);
        
        bound_tcr = interp1(time, bound_tcr, final_time, 'previous');
        bound_np = interp1(time, bound_np, final_time, 'previous');
        phos_tcr = interp1(time, phos_tcr, final_time, 'previous');
        
        homo_bound_tcr{nsims} = [bound_tcr; final_time];
        homo_bound_np{nsims} = [bound_np; final_time];
        homo_phos_tcr{nsims} = [phos_tcr; final_time];
    end
    save(ufile,'homo_bound_tcr','homo_bound_np','homo_phos_tcr');
%}
end