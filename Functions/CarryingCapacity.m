function [cnp] = CarryingCapacity(np_params, tcr_params, plt)

    % NP
    vh = np_params(1);
    rNP = np_params(2);
    np_rho = np_params(3);
    np_num = 5000;

    % TCR Params
    rSurf = tcr_params(1);
    num_clusters = tcr_params(2);
    cluster_radius = tcr_params(3);
    tcr_per_cluster = tcr_params(4);
    num_tcr = tcr_params(5);
    
    rTCR = 5;
    
    tcr_pos = gen_clusters(rSurf, cluster_radius, num_clusters, tcr_per_cluster, num_tcr);  
    
    nsims = 500;
    [cnp np_pos] = capacityNP(tcr_pos, rNP, nsims);
    
    if plt == 1
        figure()
        angles = linspace(0,2*pi,500);
        plot(rSurf*cos(angles), rSurf*sin(angles),'k-','Linewidth',2)
        hold on
        for i = 1:length(tcr_pos)
            plot(rTCR*cos(angles)+tcr_pos(1,i),rTCR*sin(angles)+tcr_pos(2,i),'b-')
        end
        for i = 1:size(np_pos,2)
            plot(rNP*cos(angles)+np_pos(1,i), rNP*sin(angles)+np_pos(2,i),'r-')
        end
    end
    
end

%% Generating random positions of points on disk
function pos = generate_pos(rSurf, dmin, num_points);
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
        %lambda = tcr_per_cluster;
        temp_pos = generate_pos(cluster_radius, rTCR, lambda);
        temp_pos = temp_pos + nc_pos(:,i);
        tcr_pos = [tcr_pos, temp_pos];
    end
    
    
    if num_clusters == 0
        tcr_pos = generate_pos(rSurf, rTCR, num_tcr);
    end
    
    %{
    free_tcr = num_tcr -length(tcr_pos);
    disp(free_tcr);
    
    if free_tcr < 0
        tcr_pos = tcr_pos(:,1:num_tcr);
    else
        temp_pos = generate_pos(rSurf, rTCR, free_tcr);
        tcr_pos = [tcr_pos, temp_pos];
    end
    %}
end
%% Compute the number of covered TCRs for each NP
function nt = covered_tcrs(tcr_pos, np_pos, rNP)
    nt = zeros(1,size(np_pos,2));
    for i = 1:size(np_pos,2)
        d = dist(np_pos(:,i)',tcr_pos);
        nt(i) = length(d(d<rNP));
    end
end

%% Estimate NP capacity
function [cnp, np_pos] = capacityNP(tcr_pos, rNP, nsims)
i=1;
np_pos=[tcr_pos(:,1)];
for i = 2:length(tcr_pos)
    
    rho = rand(1) + rand(1);
    rho(rho>1)=2-rho(rho>1);
    theta = 2*pi*rand(1);
    pos = [rNP*rho.*cos(theta); rNP*rho.*sin(theta)];
    
    np_pos = [np_pos, tcr_pos(:,i)+pos];
    
    d = dist(np_pos(:,end)',np_pos(:,1:end-1));
    min_d = min(d);
    
    j=1;
    while min_d < 2*rNP && j < nsims
        np_pos(:,end) = [];
        
        rho = rand(1) + rand(1);
        rho(rho>1)=2-rho(rho>1);
        theta = 2*pi*rand(1);
        pos = [rNP*rho.*cos(theta); rNP*rho.*sin(theta)];
    
        np_pos = [np_pos, tcr_pos(:,i)+pos];
        
        d = dist(np_pos(:,end)',np_pos(:,1:end-1));
        min_d = min(d);
        
        j = j + 1;
    end
    
    if min_d < 2*rNP && j == nsims
        np_pos(:,end) = [];
    end
    
end

cnp = length(np_pos);

end