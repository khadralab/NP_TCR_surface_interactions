function [bound_tcr, bound_np, phos_tcr, time] = Gillespie_KPR(tcr_params, np_params, sim_params)
    
    rSurf = tcr_params(1);
    num_clusters = tcr_params(2);
    cluster_radius = tcr_params(3);
    tcr_per_cluster = tcr_params(4);
    num_tcr = tcr_params(5);
    rTCR = tcr_params(6);
    
    rNP = np_params(1);
    vh = np_params(2);
    k0 = np_params(3);
    kon = np_params(4);
    koff = np_params(5);
    T0 = np_params(6);
    rho = np_params(7);
    
    kp=0.3;
    maxP = 3;

    total_t = sim_params(1);
    
    tcr_pos = gen_clusters(rSurf, cluster_radius, num_clusters, tcr_per_cluster, num_tcr);
    sample_pos = generate_NPs(rSurf, 0, total_t, []);
    sample_nt = count_covered_tcrs(tcr_pos, sample_pos, rNP);
    
    tcr_states = zeros(size(tcr_pos,2),1);
    
    t = 0;
    time = [0];
    np_pos = [];
    nt = [];
    bt = [];

    bound_tcr = [0];
    bound_np = [0];
    phos_tcr = [0];
    j=1;

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
        
        last_state = tcr_states(:,end);
        unphosphorylated = length(last_state(last_state>0 & last_state<maxP));
        ph = kp*unphosphorylated;
        

        a2 = sum(f);                % Cross-link binding propensity
        a3 = sum(b);                % TCR unbnding propensity
        a4 = sum(ph);               % TCR Phosphorylation propensity

        a0 = a1 + a2 +a3+a4;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Time of next reaction
        tau = -1/a0 * log(r1);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Update state variables
        t = t+tau;
        time = [time,t];
        tcr_states = [tcr_states, tcr_states(:,end)];
        

        % Create new sample NP positions if necessary
        if j > length(sample_pos)
            sample_pos = generate_NPs(rSurf, 0, total_t, []);
            sample_nt = count_covered_tcrs(tcr_pos, sample_pos, rNP);
            j=1;
        end

        % Probability of new NP binding
        pbind = k0*vh*sample_nt(j) / (1 + k0*vh*sample_nt(j));

        % Distance to nearest NP
        if length(np_pos) > 1
            d = dist(np_pos', sample_pos(:,j));
            min_d = min(d);
        else
            min_d = 3*rNP;
        end

        % New NP Adsorption Event
        if (r2 <= a1/a0) && (r3 < pbind) && (min_d > 2*rNP)
            np_pos = [np_pos, sample_pos(:,j)];
            nt = [nt, sample_nt(j)];
            bt = [bt, 1];
            
            tcr_states(:,end) = tcr_states(:,end-1) + adsorption(tcr_pos, sample_pos(:,j), rNP);
            j=j+1;

        % New Cross-Linking Binding    
        elseif (r2 > a1/a0) && (r2 <= (a1+a2)/a0)
            cf = cumsum(f) ./ a2;
            cf = r3 < cf;
            ind = find(cf==1,1,'first');
            bt(ind) = bt(ind)+1;
            
            tcr_states(:,end) = crosslink(tcr_pos, np_pos(:,ind), rNP, tcr_states(:,end-1));

        % TCR Unbinding Event    
        elseif (r2 > (a1+a2)/a0) && (r2 <= (a1+a2+a3)/a0)
            cb = cumsum(b) ./ a3;
            cb = r3 < cb;
            ind = find(cb == 1,1,'first');
            bt(ind) = bt(ind)-1;
            
            tcr_states(:,end) = unlink(tcr_pos, np_pos(:,ind), rNP, tcr_states(:,end-1));

            % NP Unbinding Event
            if bt(ind)==0
                tcr_states(:,end) = desorption(tcr_pos, np_pos(:,ind), rNP, tcr_states(:,end-1));
                np_pos(:,ind) = [];
                bt(ind) = [];
                nt(ind) = [];
            end
            
        % TCR Phosphorylation event
        elseif (r2 > (a1+a2+a3)/a0)
            
            tcr_states(:,end) = phosphorylation(tcr_states(:,end-1), maxP);
            
        end
        j=j+1;
        
        current_state = tcr_states(:,end);
        
        bound_tcr = [bound_tcr, sum(bt)];
        bound_np = [bound_np, size(np_pos,2)];
        phos_tcr = [phos_tcr, sum(current_state==maxP)];
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot Time Series
    %{
    g=figure();
    g.Visible = 'on';
    subplot(2,2,1)
    plot(time,bound_tcr,'b-','DisplayName',' Bound TCRs')
    xlabel('Time')
    ylabel('Bound TCRs')
    
    subplot(2,2,2)
    plot(time,bound_np,'r-','DisplayName','Bound NPs')
    xlabel('Time')
    ylabel('Bound NPs')
    
    subplot(2,2,3)
    histogram(bound_tcr,'Normalization','pdf')
    xlabel('Bound TCRs')
    ylabel('PDF')

    subplot(2,2,4)
    histogram(bound_np,'Normalization','pdf')
    xlabel('Bound NPs')
    ylabel('PDF')
    
    sgtitle([fname(1:4),': Koff = ',num2str(koff),', v = ',num2str(vh),', r = ',num2str(rNP)]);
    
    savefig(g,fname);
    %}
    
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
        d = dist(pos(:,i)',pos(:,1:i-1));
        min_d = min(d);
        nsims = 0;
        while min_d < 2*dmin && nsims < 500
            rho = rand(1,1)+rand(1,1);
            rho(rho>1) = 2-rho(rho>1);
            th = 2*pi*rand(1,1);
            pos(:,i) = [rSurf*rho*cos(th); rSurf*rho*sin(th)];
        
            d = dist(pos(:,i)',pos(:,1:i-1));
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
    
    d = dist(pos', new_pos);
    
    if d < 2*dmin
        new_pos = [5000;5000];
    end
    
    pos = [pos, new_pos];

end

%% Count the number of covered TCRs for each NP
function nt = count_covered_tcrs(tcr_pos, np_pos, rNP)
    nt = zeros(1,size(np_pos,2));
    for i = 1:size(np_pos,2)
        d = dist(np_pos(:,i)',tcr_pos);
        nt(i) = length(d(d<rNP));
    end
end

%% Return indices of TCRs covered by NP Adsorption event
function nt = adsorption(tcr_pos, np_pos, rNP)
    nt = zeros(size(tcr_pos,2),1);
    d = dist(np_pos', tcr_pos);
    
    nt(d<rNP) = -1;
    nt(find(nt == -1,1,'first')) = 1;
    
end

%% State Update: Crosslinking event
function tcr_states = crosslink(tcr_pos, np_pos, rNP, tcr_states)
    d = dist(np_pos', tcr_pos)';
    ind = find(tcr_states<0 & d<rNP);
    
    if length(ind)>1
        ind = randsample(ind,1);
    end
    
    tcr_states(ind) = 1;
    
end

%% State Update: Unlinking event
function tcr_states = unlink(tcr_pos, np_pos, rNP, tcr_states)
    d = dist(np_pos', tcr_pos)';
    ind = find(tcr_states>0 & d<rNP);
    
    if length(ind)>1
        ind = randsample(ind,1);
    end

    tcr_states(ind) = -1;
    
end

%% State Update: NP Desorption
function tcr_states = desorption(tcr_pos, np_pos, rNP, tcr_states)
    d = dist(np_pos', tcr_pos)';
    
    tcr_states(d<rNP) = 0;
    
end

%% State Update: TCR Phosphorylation
function tcr_states = phosphorylation(tcr_states,maxP)
    ind = find(tcr_states>0 & tcr_states<maxP);
    
    if length(ind) > 1
        ind = randsample(ind,1);
    end
    
    tcr_states(ind) = tcr_states(ind)+1;
end