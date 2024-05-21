function y = EMaxPlots(ax, kd_vals, rho_vals, surface_type, bound_type, colors, file)
    x_vals = linspace(min(rho_vals), max(rho_vals), 1000);
    rho_vals = 10.^rho_vals;
    x_vals = 10.^x_vals;
    
    eMax = []; y_vals = [];
    
    if strcmp(bound_type,'bound')
        k=1;
    elseif strcmp(bound_type,'phos')
        k=2;
    else
        k=3;
    end
    
    i=1;
    
    for koff = kd_vals ./ 10
        [mean_data, var_data, output_fit] = load_DR(koff, rho_vals, surface_type, file);
        fit_data = output_fit{k};
        eMax = [eMax,fit_data.eMax];
        y_vals = [y_vals, fit_data(x_vals)];
    end
    
    x = log10(kd_vals); y = eMax';
    X = [ones(length(x),1),x'];

    b = X \ y;
    
    hold(ax,'on');
    scatter(ax, kd_vals, eMax, [], colors, 'filled');
    %plot(ax, kd_vals, 10^b(1)*kd_vals.^b(2));
    ylabel(ax,'EMax');
    xlabel(ax,'$K_D$');
    grid(ax,'on');
    
end