function [bound_tcr, bound_np, phos_tcr, tcr_states, time] = Gillespie_v3(tcr_params, np_params, sim_params)
    
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

    total_t = sim_params(1);
    
    num_samples = floor(100);
    
    % Generate TCR positions and pre-define NP positions to sample from.
    % Determine covered tcrs (nt) for each NP pos.
    tcr_pos = gen_clusters(rSurf, cluster_radius, num_clusters, tcr_per_cluster, num_tcr);
    sample_pos = generate_NPs(rSurf, 0, num_samples, []);
    
    tcr_pos = single(tcr_pos); sample_pos = single(sample_pos);
    
    sample_nt = count_covered_tcrs(tcr_pos, sample_pos, rNP);
    dist_tcr_to_np = mdist(sample_pos, tcr_pos);          % Distance of each possible NP to each TCR: size=(num_NP, num_tcr)
    
    % Initialize state variables
    tcr_states = zeros(1,size(tcr_pos,2));
    
    t = 0;                  % Current simulation time
    time = [0];             % Array of reaction times
    np_pos = 5000*ones(2,300);        % Array of bound NP positions
    nt = zeros(1,300);                % Current # TCRs covered by each bound NP
    bt = zeros(1,300);                % Current # TCRs bound by each bound NP

    id_range = [1:300];             % Range of possible simultaneously bound NPs; Max = 300.
    tcr_np_id = zeros(1, num_tcr);  % For each TCR identify the NP index which is covering that TCR.
    
    bound_tcr = [0];        % Array of bound TCRs over time
    bound_np = [0];         % Array of bound NPs over time
    phos_tcr = [0];         % Array of Phos TCRs over time
    dist_tcr_to_bound_np = mdist(np_pos, tcr_pos);
    
    j=1;                    % Index of NP to sample from

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Gillespie Algorithm

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
        %tcr_states = [tcr_states, tcr_states(:,end)];

        % Create new NP positions to sample from, if necessary.
        if j > length(sample_pos)
            sample_pos = generate_NPs(rSurf, 0, num_samples, []);
            sample_pos = single(sample_pos);
            sample_nt = count_covered_tcrs(tcr_pos, sample_pos, rNP);
            dist_tcr_to_np = mdist(sample_pos, tcr_pos);
            
            j=1;
        end

        % Probability of new NP binding
        pbind = k0*vh*sample_nt(j) / (1 + k0*vh*sample_nt(j));

        % Distance to nearest NP
        if length(np_pos) > 1
            d = mdist(np_pos, sample_pos(:,j));
            min_d = min(d);
        else
            min_d = 3*rNP;
        end

        % New NP Adsorption Event
        if (r2 <= a1/a0) && (min_d > 2*rNP) && (r3 < pbind)
            next_np_id = ismember(id_range, tcr_np_id);
            next_np_id = id_range(~next_np_id);
            next_np_id = next_np_id(1);
            
            np_pos(:,next_np_id) = sample_pos(:,j);
            
            np_pos = single(np_pos); tcr_pos = single(tcr_pos);             % Convert positions to single precision to save memory.
            dist_tcr_to_bound_np(next_np_id,:) = mdist(np_pos(:,next_np_id), tcr_pos);
            
            tcr_np_id(dist_tcr_to_bound_np(next_np_id,:) < rNP) = next_np_id;
            
            tcr_states(tcr_np_id == next_np_id) = -1;
            
            check = find(tcr_np_id == next_np_id);
            
            if isempty(check)
                disp(['Next NP ID:', num2str(next_np_id)]);
                disp(['Total Bound NPs:', num2str(max(tcr_np_id))]);
                error('Error: Cannot sample from empty array!');
            else
                tcr_indx = datasample(check, 1,'Replace',false);
                tcr_states(tcr_indx) = 1;
            end
            
            nt(next_np_id) = sum(tcr_np_id == next_np_id);
            bt(next_np_id) = 1;
            
            % Check sum and raise error
            if sum(bt) ~= sum(tcr_states>0)
                error('Bound TCRs do not agree')
            end
            
        % New Cross-Linking Binding    
        elseif (r2 > a1/a0) && (r2 <= (a1+a2)/a0)
            cf = cumsum(f) ./ a2;
            cf = r3 < cf;
            np_indx = find(cf==1,1,'first');
            bt(np_indx) = bt(np_indx)+1;
            
            condition = find(tcr_states < 0 & tcr_np_id == np_indx);
            tcr_indx = datasample(condition, 1,'Replace', false);
            
            tcr_states(tcr_indx) = 1;
            
            % Check sum and raise error
            if sum(bt) ~= sum(tcr_states>0)
                error('Bound TCRs do not agree')
            end
            
            j=j-1;

        % TCR Unbinding Event    
        elseif (r2 > (a1+a2)/a0) && (r2 <= (a1+a2+a3)/a0)
            cb = cumsum(b) ./ a3;
            cb = r3 < cb;
            np_indx = find(cb == 1,1,'first');
            bt(np_indx) = bt(np_indx)-1;
            
            tcr_indx = datasample(find(tcr_states > 0 & tcr_np_id == np_indx), 1, 'Replace', false);
            tcr_states(tcr_indx) = -1;

            % NP Unbinding Event
            if bt(np_indx)==0
                tcr_states(tcr_np_id == np_indx) = 0;
                tcr_np_id(tcr_np_id == np_indx) = 0;
                
                np_pos(:,np_indx) = [5000; 5000];
                dist_tcr_to_bound_np(np_indx,:) = 5000*ones(1,300);

                bt(np_indx) = 0;
                nt(np_indx) = 0;
            end
            
            % Check sum and raise error
            if sum(bt) ~= sum(tcr_states>0)
                error('Bound TCRs do not agree')
            end
            j=j-1;
            
        % TCR Phosphorylation event
        elseif (r2 > (a1+a2+a3)/a0)
            
            tcr_indx = find(tcr_states>0 & tcr_states<maxP);
    
            if ~isempty(tcr_indx)
                tcr_indx = datasample(tcr_indx,1,'Replace',false);
                tcr_states(tcr_indx) = tcr_states(tcr_indx)+1;
            end
            
            % Check sum and raise error
            if sum(bt) ~= sum(tcr_states>0)
                error('Bound TCRs do not agree')
            end
            
            j=j-1;
            
        end
        
        j=j+1;
        
        bound_tcr = [bound_tcr, sum(bt)];
        bound_np = [bound_np, max(tcr_np_id)];
        phos_tcr = [phos_tcr, sum(tcr_states==maxP)];
    end
    
end

%% Local Functions

%% Generating random positions of points on disk
function pos = generate_pos(rSurf, dmin, num_points)
    rho = rand(1,num_points) + rand(1,num_points);
    rho(rho>1)=2-rho(rho>1);
    theta = 2*pi*rand(1,num_points);

    pos = [rSurf*rho.*cos(theta); rSurf*rho.*sin(theta)];

    i=2;

    while i<=length(pos)
        d = mdist(pos(:,i),pos(:,1:i-1));
        min_d = min(d);
        nsims = 0;
        while min_d < 2*dmin && nsims < 500
            rho = rand(1,1)+rand(1,1);
            rho(rho>1) = 2-rho(rho>1);
            th = 2*pi*rand(1,1);
            pos(:,i) = [rSurf*rho*cos(th); rSurf*rho*sin(th)];
        
            d = mdist(pos(:,i),pos(:,1:i-1));
            min_d = min(d);
            nsims = nsims+1;
        end
        
        if nsims == 500
            pos(:,i)=[];
            i=i-1;
        end
        i=i+1;
    end
end
%
%
%
%% Generate TCR Nanoclusters (Alternative Gaussian Filtering)
function tcr_pos = gen_clusters(rSurf, cluster_radius, num_clusters, tcr_per_cluster, num_tcr)
    rTCR = 5;
    nc_pos = generate_pos(rSurf-cluster_radius, 2*cluster_radius, num_clusters);
    tcr_pos = [];
    for i = 1:num_clusters
        lambda = poissrnd(tcr_per_cluster);
        lambda = tcr_per_cluster;
        temp_pos = generate_pos(cluster_radius, rTCR, lambda);
        temp_pos = temp_pos + nc_pos(:,i);
        tcr_pos = [tcr_pos, temp_pos];
    end
    
    free_tcr = num_tcr -length(tcr_pos);
    
    if free_tcr < 0
        tcr_pos = tcr_pos(:,1:num_tcr);
    else
        temp_pos = generate_pos(rSurf, rTCR, free_tcr);
        tcr_pos = [tcr_pos, temp_pos];
    end
end

%% Generating random NP positions
function pos = generate_NPs(rSurf, dmin, num_points, pos)
    rho = rand(1,num_points) + rand(1,num_points);
    rho(rho>1)=2-rho(rho>1);
    theta = 2*pi*rand(1,num_points);

    new_pos = [rSurf*rho.*cos(theta); rSurf*rho.*sin(theta)];
    
    %{
    d = dist(pos', new_pos);
    
    if d < 2*dmin
        new_pos = [5000;5000];
    end
    %}
    pos = [pos, new_pos];

end

%% Count the number of covered TCRs for each NP
function nt = count_covered_tcrs(tcr_pos, np_pos, rNP)
    d = mdist(np_pos, tcr_pos);
    nt = sum(d<rNP,2);
end

%% Custom Distance Function

function D = mdist(X,Y)
    D = bsxfun(@plus,dot(X,X,1)',dot(Y,Y,1))-2*(X'*Y);
    D = sqrt(D);
end