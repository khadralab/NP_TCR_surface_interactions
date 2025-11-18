function [bound_tcr, bound_np, phos_tcr, tcr_states, time] = ClusterTau(tcr_params, np_params, sim_params)
    
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
    tau = sim_params(2);
    window = 1000000;
    thresh = 1;
    
    num_samples = floor(100);
    
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
    
    j=1;                                    % Index of NP to sample from

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Tau Leaping Algorithm
    convergence = false;
    
    while t < total_t
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Calculate Propensities
        a1 = rho*T0;                % NP adsorption propensity

        f = kon*(vh-bt).*(nt-bt);
        b = koff*bt;

        unphosphorylated =  sum((tcr_states>0) & (tcr_states<maxP), 'all');
        ph = kp*unphosphorylated;
        

        a2 = sum(f);                % Cross-link binding propensity
        a3 = sum(b);                % TCR unbinding propensity
        a4 = sum(ph);               % TCR Phosphorylation propensity

        a0 = a1+a2+a3+a4;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Time of next reaction
        tau = sim_params(2);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Frequency of each reaction within time-step tau.
        
        e1 = poissrnd(a1*tau);      % Number of new NP arrivals
        e2 = poissrnd(a2*tau);      % Number of cross-linking reactions
        e3 = poissrnd(a3*tau);      % Number of unbinding events
        e4 = poissrnd(a4*tau);      % Number of phosphorylation events
        
        E = [e2, e3, e4];
        
        while any(E > bound_np(end))           % Condition to avoid large time steps with too many reactions.
            tau = tau * bound_np(end) / max(E);
            
            e1 = poissrnd(a1*tau);      % Number of new NP arrivals
            e2 = poissrnd(a2*tau);      % Number of cross-linking reactions
            e3 = poissrnd(a3*tau);      % Number of unbinding events
            e4 = poissrnd(a4*tau);      % Number of phosphorylation events

            E = [e2, e3, e4];
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Update state variables
        t = t+tau;
        time = [time,t];
        %tcr_states = [tcr_states, tcr_states(:,end)];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % NP adsorption        
        if e1 ~= 0
            [tcr_states, np_pos] = updateBoundNPs(tcr_states, tcr_pos, np_pos, e1, rSurf, rNP, k0, vh);
        end    
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % TCR cross-Linking
        
        if e2 ~= 0
            [tcr_states] = crosslinking(tcr_states,f,e2,vh,kon);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % TCR un-linking
        
        if e3 ~= 0
            [tcr_states, bt, nt] = unlinking(tcr_states,b,e3,koff);
        end
        
        removeNP = find((bt==0) & (nt~=0));
        
        if ~isempty(removeNP)
            tcr_states(:,removeNP) = zeros(num_tcr,length(removeNP));
            np_pos(:,removeNP) = 5000*ones(2,length(removeNP));
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % TCR Phosphorylation
        
        if e4 ~= 0
            [tcr_states] = phosphorylation(tcr_states, e4, maxP);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Count Bound Receptors and NPs
        
        bt = sum(tcr_states > 0, 1);
        nt = sum(tcr_states ~= 0, 1);
        
        bound_tcr = [bound_tcr, sum(tcr_states>0,'all')];
        bound_columns = sum(abs(tcr_states),1);
        bound_np = [bound_np, sum(bound_columns > 0)];
        phos_tcr = [phos_tcr, sum(tcr_states == 3,'all')];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Check convergence of state variables
        
        % Identify slowest kinetic variable
        %{
        slow_var = min([kon, koff, kp, k0, rho*T0]);
        window = floor(1/slow_var);
        
        conv1 = checkconvergence(bound_tcr, window, tau);
        conv2 = checkconvergence(bound_np, window, tau);
        conv3 = checkconvergence(phos_tcr, window, tau);
        
        if conv1 == true && conv2 == true && conv3 == true
            convergence = true;
        end
        %}
    end
end
