function [tcr_pos, np_pos, cnp] = NPcapacity(rSurf, cluster_radius, num_clusters, tcr_per_cluster, num_tcr, rNP)
    
    tcr_pos = generate_pos(cluster_radius, rTCR, num_tcr);
    
    nsims= 0;
    
    while free_tcr > 0 && nsims <500
        rho = rand(1,1) + rand(1,1);
        rho(rho>1)=2-rho(rho>1);
        theta = 2*pi*rand(1,1);
        
        new_np_pos = 

        pos = [rSurf*rho.*cos(theta); rSurf*rho.*sin(theta)];

    
        dist_tcr_np = dist(new_np_pos', tcr_pos);
    
        nt = sum(dist_np_tcr<rNP,2)';
        
        nsims= nsims+1;
    
end

