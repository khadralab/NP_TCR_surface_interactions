function [cluster_states] = MixtureResponseFunction(rho_index, mixRatio)
    
    % NP and Kinetic Params
    % Strong ligand: 1-10s --> koff = 0.1-1
    % Weak ligand: 0.1-3s --> koff = 0.3-10
    
    % NP 1 Params
    rNP1 = 10;
    vh1 = 2;
    k01 = 0.1;
    kon1 = 0.1;
    koff1 = 0.01;
    T01 = 1;

    % NP 2 Params
    rNP2 = 10;
    vh2 = 2;
    k02 = 0.1;
    kon2 = 0.1;
    koff2 = 0.1;
    T02 = 1;
    
    
    % Sim Params
    computecanada = false;
    total_t = 1e6;
    tau = 1;
    sim_params = [total_t, tau];
    total_sims = 10;
    
    % NP concentrations
    rho_vals = linspace(-2,5,15);
    rho_vals = 10.^rho_vals;
    
    rho_total = rho_vals(rho_index);
    
    rho1 = rho_total * mixRatio;
    rho2 = rho_total * (1-mixRatio);
    
    np1_params = [rNP1, vh1, k01, kon1, koff1, T01, rho1];
    np2_params = [rNP2, vh2, k02, kon2, koff2, T02, rho2];

    % Cluster Surface
    rSurf = 1000;
    num_clusters = 15;
    cluster_radius = 50;
    tcr_per_cluster = 20;
    num_tcr = 300;
    rTCR = 5;
    tcr_params = [rSurf, num_clusters, cluster_radius, tcr_per_cluster, num_tcr, rTCR];

    cfile = ['MixtureResults/temp/Cluster_rho',num2str(rho_total),'-',num2str(mixRatio*100),'.mat'];
    %ufile = ['MixtureResults/Uniform_rho',num2str(rho),'-',num2str(mixRatio*100),'.mat'];
    
    if computecanada == true
         % Create a "local" cluster object
        local_cluster = parcluster('local');

        % Modify the JobStorageLocation to $SLURM_TMPDIR
        local_cluster.JobStorageLocation = getenv('SLURM_TMPDIR');

        % Start the parallel pool
        parpool(local_cluster)
        fprintf("Successfully initiated parpool \n")
        
        cfile = ['MixtureResults/Cluster_rho',num2str(rho),'-',num2str(mixRatio*100),'.mat'];
    end
    
    parfor nsims = [1:total_sims]
        [bound_tcr, bound_np1, bound_np2, phos_tcr, cluster_states, time] = TauMixtures(tcr_params, np1_params, np2_params, sim_params);
        [time,j,~] = unique(time);
        total_t = floor(time(end))+1;
        %final_time = linspace(0,total_t,total_t*10);

        bound_tcr = bound_tcr(j); bound_np1 = bound_np1(j); bound_np2 = bound_np2(j); phos_tcr = phos_tcr(j);
        
        %bound_tcr = interp1(time, bound_tcr, final_time, 'previous');
        %bound_np = interp1(time, bound_np, final_time, 'previous');
        %phos_tcr = interp1(time, phos_tcr, final_time, 'previous');
        
        cluster_bound_tcr{nsims} = [bound_tcr; time];
        cluster_bound_np1{nsims} = [bound_np1; time];
        cluster_bound_np2{nsims} = [bound_np2; time];
        cluster_phos_tcr{nsims} = [phos_tcr; time];
    end
    
    save(cfile,'cluster_bound_tcr','cluster_bound_np1','cluster_bound_np2','cluster_phos_tcr');
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