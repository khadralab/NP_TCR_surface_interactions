 function tcr_pos = gen_clusters(rSurf, cluster_radius, num_clusters, tcr_per_cluster, num_tcr)
    rTCR = 5;
    nc_pos = generate_pos(rSurf-cluster_radius, 2*cluster_radius, num_clusters);
    tcr_pos = [];
    
    if tcr_per_cluster < 1
         lambda = poissrnd(tcr_per_cluster,1,num_clusters);
    else
         lambda = tcr_per_cluster*ones(1, num_clusters);
    end
        
    
    for i = 1:num_clusters
        temp_pos = generate_pos(cluster_radius, rTCR, lambda(i));
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