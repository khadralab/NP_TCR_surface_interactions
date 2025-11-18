function [time_cluster, time_uniform, cluster_states, uniform_states] = DoseResponseUniform(rho_index)
    
    % NP and Kinetic Params
    % Strong ligand: 1-10s --> koff = 0.1-1
    % Weak ligand: 0.1-3s --> koff = 0.3-10
    rNP = 20;
    vh = 5;
    k0 = 0.1;
    kon = 0.1;
    koff = 0.1;
    T0 = 1;

    % Simulation Params
    computecanada = true;
    
    if rho_index > 70
        total_t = 1e5;
        tau = 500;
        sim_params = [total_t, tau];
    else
        total_t = 1e6;
        tau = 500;
        sim_params = [total_t, tau];
    end
    total_sims = 30;
    
    
    % NP concentrations
    rho_vals = linspace(-2,5,15);
    rho_vals = 10.^rho_vals;
    rho = rho_vals(rho_index);
    
    np_params = [rNP, vh, k0, kon, koff, T0, rho];

    fprintf(['rho = ',num2str(rho),' \n']);
    
    % Homogeneous Surface
    rSurf = 1000;
    num_clusters = 0;
    cluster_radius = 50;
    tcr_per_cluster = 20;
    num_tcr = 300;
    rTCR = 5;
    tcr_params = [rSurf, num_clusters, cluster_radius, tcr_per_cluster, num_tcr, rTCR];

    ufile = ['LongSims/temp/koff',num2str(koff*1000),'/Uniform_rho',num2str(floor(rho*100))];
    
    if computecanada == true
         % Create a "local" cluster object
        local_cluster = parcluster('local');

        % Modify the JobStorageLocation to $SLURM_TMPDIR
        local_cluster.JobStorageLocation = getenv('SLURM_TMPDIR');

        % Start the parallel pool
        parpool(local_cluster)
        fprintf("Successfully initiated parpool \n")
        
        ufile = ['koff',num2str(koff*100),'/Uniform_rho',num2str(floor(rho*100))];
    end
    
    %save(ufile,'homo_bound_tcr','homo_bound_np','homo_phos_tcr');
    
    parfor nsims = 1:total_sims
        [bound_tcr, bound_np, phos_tcr, cluster_states, time] = TauLeaping(tcr_params, np_params, sim_params);
        [time,j,~] = unique(time);
        total_t = floor(time(end))+1;

        bound_tcr = bound_tcr(j); bound_np = bound_np(j); phos_tcr = phos_tcr(j);
        
        homo_bound_tcr{nsims} = [bound_tcr; time];
        homo_bound_np{nsims} = [bound_np; time];
        homo_phos_tcr{nsims} = [phos_tcr; time];
    end
    
    if ~exist(['LongSims/temp/koff',num2str(koff*100),'/'])
        mkdir(['LongSims/temp/koff',num2str(koff*100),'/']);
    end
    save(ufile,'homo_bound_tcr','homo_bound_np','homo_phos_tcr');

end
