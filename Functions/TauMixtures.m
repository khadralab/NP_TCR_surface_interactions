function [bound_tcr, bound_np1, bound_np2, phos_tcr, tcr_states, time] = TauMixtures(tcr_params, np1_params, np2_params, sim_params)
    
    addpath('Functions')

    % TCR params
    rSurf = tcr_params(1);
    num_clusters = tcr_params(2);
    cluster_radius = tcr_params(3);
    tcr_per_cluster = tcr_params(4);
    num_tcr = tcr_params(5);
    
    % NP 1 params
    rNP1 = np1_params(1);
    vh1 = np1_params(2);
    k01 = np1_params(3);
    kon1 = np1_params(4);
    koff1 = np1_params(5);
    T01 = np1_params(6);
    rho1 = np1_params(7);
    
    % NP 1 params
    rNP2 = np2_params(1);
    vh2 = np2_params(2);
    k02 = np2_params(3);
    kon2 = np2_params(4);
    koff2 = np2_params(5);
    T02 = np2_params(6);
    rho2 = np2_params(7);
    
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
    np_type = zeros(1, num_tcr);            % 0-No NPs; 1-First NP; 2-Second NP
    
    %[tcr_states, np_pos, bound_tcr, covered_tcr] = initialConditions(tcr_states, tcr_pos, np_pos, rSurf, rNP, vh);
            
    t = 0;                              % Current simulation time
    time = t;                           % Array of reaction times
    
    nt = zeros(num_tcr,1) ;        % Covered TCRs for each NP type 1
    bt = nt;                      % Bound TCRs for each NP type 1
    
    bound_tcr = sum(bt);                          % Array of bound TCRs over time
    bound_np1 = sum(np_pos(1,np_type==1) ~= 5000);                 % Array of bound NPs over time
    bound_np2 = sum(np_pos(1,np_type==2) ~= 5000);
    phos_tcr = sum(tcr_states == 3, 'all');              % Array of Phos TCRs over time
    
    j=1;                                    % Index of NP to sample from

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Tau Leaping Algorithm
    convergence = false;
    
    while t < total_t
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Calculate Propensities
        a1 = rho1*T01;
        a7 = rho2*T02;                % NP adsorption propensity

        f1 = kon1*(vh1-bt(np_type==1)).*(nt(np_type==1)-bt(np_type==1));
        b1 = koff1*bt(np_type==1);
        
        f2 = kon2*(vh2-bt(np_type==2)).*(nt(np_type==2)-bt(np_type==2));
        b2 = koff2*bt(np_type==2);

        unphosphorylated =  sum((tcr_states>0) & (tcr_states<maxP), 'all');
        ph = kp*unphosphorylated;
        

        a2 = sum(f1);                % Cross-link binding propensity
        a3 = sum(b1);                % TCR unbinding propensity
        a4 = sum(ph);                % TCR Phosphorylation propensity
        
        a5 = sum(f2);
        a6 = sum(b2);

        a0 = a1+a2+a3+a4+a5+a6+a7;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Time of next reaction
        tau = sim_params(2);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Frequency of each reaction within time-step tau.
        
        e1 = poissrnd(a1*tau);      % Number of new NP arrivals
        e2 = poissrnd(a2*tau);      % Number of cross-linking reactions
        e3 = poissrnd(a3*tau);      % Number of unbinding events
        e4 = poissrnd(a4*tau);      % Number of phosphorylation events
        
        e5 = poissrnd(a5*tau);      % Number of cross-linking for NP species 2
        e6 = poissrnd(a6*tau);      % Number of unbinding for NP species 2
        e7 = poissrnd(a7*tau);
        
        E = [e2, e3, e4, e5, e6];
        
        while any(E > (bound_np1(end)+bound_np2(end)))           % Condition to avoid large time steps with too many reactions.
            tau = tau * (bound_np1(end)+bound_np2(end)) / max(E);
            
            e1 = poissrnd(a1*tau);      % Number of new NP arrivals
            e2 = poissrnd(a2*tau);      % Number of cross-linking reactions
            e3 = poissrnd(a3*tau);      % Number of unbinding events
            e4 = poissrnd(a4*tau);      % Number of phosphorylation events
            
            e5 = poissrnd(a5*tau);      % Number of cross-linking for NP species 2
            e6 = poissrnd(a6*tau);      % Number of unbinding for NP species 2
            e7 = poissrnd(a7*tau);

            E = [e2, e3, e4, e5, e6];
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Update state variables
        t = t+tau;
        time = [time,t];
        %tcr_states = [tcr_states, tcr_states(:,end)];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % TCR cross-Linking
        %{
        if any((any(tcr_states(:, np_type==1) < 0,1)) ~= (f1~=0))
            error('Something Wrong')
        elseif any((any(tcr_states(:, np_type==2) < 0,1)) ~= (f2~=0))
            error('Something Still Wrong')
        end
        %}
        
        if e2 ~= 0
            [tcr_states(:, np_type==1)] = crosslinking(tcr_states(:, np_type==1),f1,e2,vh1,kon1);
        end
        
        if e5 ~= 0
            [tcr_states(:, np_type==2)] = crosslinking(tcr_states(:, np_type==2),f2,e5,vh2,kon2);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % TCR un-linking
        
        if e3 ~= 0
            [tcr_states(:, np_type==1), bt(np_type==1), nt(np_type==1)] = unlinking(tcr_states(:, np_type==1),b1,e3,koff1);
        end
        
        if e6 ~= 0
            [tcr_states(:, np_type==2), bt(np_type==2), nt(np_type==2)] = unlinking(tcr_states(:, np_type==2),b2,e6,koff2);
        end
        
        removeNP = find((bt==0) & (nt~=0));
        
        if ~isempty(removeNP)
            tcr_states(:,removeNP) = zeros(num_tcr,length(removeNP));
            np_pos(:,removeNP) = 5000*ones(2,length(removeNP));
            np_type(removeNP) = 0;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % TCR Phosphorylation
        
        if e4 ~= 0
            [tcr_states] = phosphorylation(tcr_states, e4, maxP);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % NP adsorption        
        if e1 ~= 0
            [tcr_states, np_pos, np_type] = mixtureBoundNPs(tcr_states, tcr_pos, np_pos, 1, np_type, e1, rSurf, np1_params, np2_params);
        end
        
        if e7 ~= 0
            [tcr_states, np_pos, np_type] = mixtureBoundNPs(tcr_states, tcr_pos, np_pos, 2, np_type, e7, rSurf, np1_params, np2_params);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Count Bound Receptors and NPs
        
        bt = sum(tcr_states > 0, 1);
        nt = sum(tcr_states ~= 0, 1);
        
        bound_tcr = [bound_tcr, sum(tcr_states>0,'all')];
        %bound_columns = sum(abs(tcr_states),1);
        bound_np1 = [bound_np1, sum(np_type==1)];
        bound_np2 = [bound_np2, sum(np_type==2)];
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
