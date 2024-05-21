function [alpha] = ClusterDiscrimination(ax, kd_vals, rho_vals, bound_type, col)
    
    if strcmp(bound_type,'bound')
        k=1;
    elseif strcmp(bound_type,'phos')
        k=2;
    else
        k=3;
    end

    rho_vals = (10.^rho_vals)./100;

    tcrs_per_cluster = [3,5,10];
    alpha = [];
    
    for tpc = tcrs_per_cluster
        
        if tpc == 0
            file = ['TauLeap/v5'];
            surface_type = 'uniform';
        elseif tpc == 20
            file = ['TauLeap/v5'];
            surface_type = 'clusters';
        else
            file = ['TauLeap/',num2str(tpc),'TPC'];
            surface_type = 'clusters';
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
    
    bar(ax,tcrs_per_cluster, alpha)
    xlabel(ax,'TCRs per Cluster')
    ylabel(ax,'Cooperativity')
    
end