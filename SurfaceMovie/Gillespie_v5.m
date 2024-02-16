function [bound_tcr, bound_np, phos_tcr, tcr_states, time] = Gillespie_v5(tcr_params, np_params, sim_params)

    addpath('Functions')
    
    % TCR params
    rSurf = tcr_params(1);
    num_clusters = tcr_params(2);
    cluster_radius = tcr_params(3);
    tcr_per_cluster = tcr_params(4);
    num_tcr = tcr_params(5);
    
    % NP params
    rNP = np_params(1);
    vh = np_params(2);
    k0 = np_params(3);
    kon = np_params(4);
    koff = np_params(5);
    T0 = np_params(6);
    rho = np_params(7);
    
    % Phosphorylation rates and signaling state
    kp=0.3;
    maxP = 3;
    
    % Simulation params and convergence criterion
    total_t = sim_params(1);
    window = 2000;
    thresh = 1;
    
    num_samples = floor(100);
    
    % Generate TCR positions and pre-define NP positions to sample from.
    % Determine covered tcrs (nt) for each NP pos.
    tcr_pos = gen_clusters(rSurf, cluster_radius, num_clusters, tcr_per_cluster, num_tcr);
    sample_pos = generate_NPs(rSurf, 0, num_samples, []);
    
    tcr_pos = single(tcr_pos); sample_pos = single(sample_pos);
    
    sample_nt = count_covered_tcrs(tcr_pos, sample_pos, rNP);
    %dist_tcr_to_np = mdist(sample_pos, tcr_pos);          % Distance of each possible NP to each TCR: size=(num_NP, num_tcr)
    
    % Initialize state variables
    tcr_states = zeros(1,size(tcr_pos,2));
    
    np_pos = [];
    
    bound_np = 0; init_np = floor(100*rand(1));
    
    % Generate TCR positions and pre-define NP positions to sample from.
    % Determine covered tcrs (nt) for each NP pos.
    tcr_pos = gen_clusters(rSurf, cluster_radius, num_clusters, tcr_per_cluster, num_tcr);
    
    tcr_pos = single(tcr_pos);
    
    % Initialize state variables
    tcr_states = zeros(num_tcr,num_tcr);    % Num TCR x Max Num NPs
    np_pos = 5000*ones(2,num_tcr);          % Array of bound NP positions
    np_pos = single(np_pos);
    
    [tcr_states, np_pos, bound_tcr, covered_tcr] = initialConditions(tcr_states, tcr_pos, np_pos, rSurf, rNP, vh);
            
    t = 0;                                  % Current simulation time
    time = t;                               % Array of reaction times
    
    nt = covered_tcr;                 % Covered TCRs for each NP
    bt = bound_tcr;                 % Bound TCRs for each NP
    
    bound_tcr = sum(bound_tcr);                          % Array of bound TCRs over time
    bound_np = sum(np_pos(1,:) ~= 5000);                 % Array of bound NPs over time
    phos_tcr = sum(tcr_states == 3, 'all');              % Array of Phos TCRs over time
    
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Gillespie Algorithm
    convergence = false;
    
    while t < total_t
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Generate Random Numbers
        r1 = rand(1);
        r2 = rand(1);
        r3 = rand(1);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Calculate Propensities
        a1 = rho*T0;                % NP adsorption propensity

        f = kon*(vh-bt).*(nt-bt);
        b = koff*bt;

        unphosphorylated = length(tcr_states(tcr_states>0 & tcr_states<maxP));
        ph = kp*unphosphorylated;
        

        a2 = sum(f);                % Cross-link binding propensity
        a3 = sum(b);                % TCR unbnding propensity
        a4 = sum(ph);               % TCR Phosphorylation propensity

        a0 = a1+a2+a3+a4;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Time of next reaction
        tau = -1/a0 * log(r1);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Update state variables
        t = t+tau;
        time = [time,t];
        
        % New NP Adsorption Event
        if (r2 <= a1/a0)
            [tcr_states, np_pos] = updateBoundNPs(tcr_states, tcr_pos, np_pos, 1, rSurf, rNP, k0, vh);
        
        % New Cross-Linking Binding    
        elseif (r2 > a1/a0) && (r2 <= (a1+a2)/a0)
            [tcr_states] = crosslinking(tcr_states,f,1,vh,kon);
            
        % TCR Unbinding Event    
        elseif (r2 > (a1+a2)/a0) && (r2 <= (a1+a2+a3)/a0)
            [tcr_states, bt, nt] = unlinking(tcr_states,b,1,koff);
            
            removeNP = find((bt==0) & (nt~=0));
            
            % Remove NP if no TCRs are bound to it
            if ~isempty(removeNP)
                tcr_states(:,removeNP) = zeros(num_tcr,length(removeNP));
                np_pos(:,removeNP) = 5000*ones(2,length(removeNP));
            end
            
        % TCR Phosphorylation event
        elseif (r2 > (a1+a2+a3)/a0)
            [tcr_states] = phosphorylation(tcr_states, 1, maxP);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Count Bound Receptors and NPs
        
        bt = sum(tcr_states > 0, 1);
        nt = sum(tcr_states ~= 0, 1);
        
        bound_tcr = [bound_tcr, sum(tcr_states>0,'all')];
        bound_columns = sum(abs(tcr_states),1);
        bound_np = [bound_np, sum(bound_columns > 0)];
        phos_tcr = [phos_tcr, sum(tcr_states == 3,'all')];
        
    end