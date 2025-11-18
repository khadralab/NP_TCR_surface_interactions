function [alpha] = ClusterDiscrimination(ax, kd_vals, rho_vals, valence, bound_type, col)
    
    if strcmp(bound_type,'bound')
        k=1;
    elseif strcmp(bound_type,'phos')
        k=2;
    else
        k=3;
    end

    rho_vals = 10.^rho_vals;

    tcrs_per_cluster = [0,3,5,10,20];
    alpha = [];
    
    for tpc = tcrs_per_cluster
        
        if tpc == 0
            surface_type = 'uniform';
            file = ['TauLeap/0TPC/v',num2str(valence)];
        else
            surface_type = 'clusters';
            file = ['TauLeap/',num2str(tpc),'TPC/v',num2str(valence)];
        end
        
        ec50 = [];

        for koff = kd_vals ./ 10
            [mean_data, var_data, output_fit] = load_DR(koff, rho_vals, surface_type, file);
            fit_data = output_fit{k};
            ec50 = [ec50,fit_data.ec50];
        end

        x = log10(kd_vals); y = log10(ec50)';
        X = [ones(length(x),1),x'];

        b = X \ y;

        alpha = [alpha, b(2)];
    end
    
    plot(ax,tcrs_per_cluster, alpha, ':x','Color', col, 'Linewidth',1.5)
    xlabel(ax,'TCRs per Cluster')
    ylabel(ax,'Cooperativity')
    
end