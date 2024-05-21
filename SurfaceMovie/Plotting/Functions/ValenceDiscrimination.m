function [alpha] = ValenceDiscrimination(ax, kd_vals, rho_vals, surface_type, bound_type, col)
    
    if strcmp(bound_type,'bound')
        k=1;
    elseif strcmp(bound_type,'phos')
        k=2;
    else
        k=3;
    end

    rho_vals = 10.^rho_vals;

    valences = [1,2,3,5,10];
    alpha = [];
    
    for v = valences
        file = ['TauLeap/v',num2str(v)];
        ec50 = [];

        i=1;

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
    
    plot(ax,valences, alpha, ':x', 'Color', col)
    xlabel(ax,'NP Valence')
    ylabel(ax,'Cooperativity')
    
end