function [alpha, line] = ThresholdActivation(ax, kd_vals, rho_vals, tpc, valence, bound_type, colors)
    x_vals = linspace(min(rho_vals), max(rho_vals), 1000);
    rho_vals = 10.^rho_vals;
    x_vals = 10.^x_vals;
    x2_vals = x_vals;
    
    threshold = 15;
    
    P30 = [];
    
    if tpc == 0
        surface_type = 'uniform';
        file = ['TauLeap/0TPC/v',num2str(valence)];
    elseif tpc == 20
        surface_type = 'clusters';
        file = ['TauLeap/20TPC/v',num2str(valence)];
    else
        surface_type = 'clusters';
        file = ['TauLeap/',num2str(tpc),'TPC/v',num2str(valence)];
    end
    
    i=1;
    
    for koff = kd_vals ./ 10
        [mean_data, var_data, output_fit] = load_DR(koff, rho_vals, surface_type, file,1);
        
        if strcmp(bound_type,'bound')
            data = mean_data{1};
            err = var_data{1};
            fitvals = output_fit{1};
            y_vals = fitvals(x_vals);
            xlbl = 'Nanoparticle Concentration';
        elseif strcmp(bound_type, 'rescaled')
            fitvals = output_fit{1};
            data = mean_data{1} ./ fitvals.eMax;
            y_vals = fitvals(x2_vals) ./ fitvals.eMax;
            x_vals = x2_vals ./ fitvals.ec50;
            rho_vals = np_rho ./ fitvals.ec50;
            xlbl = 'Nanoparticle Concentration';
        elseif strcmp(bound_type,'phos')
            data = mean_data{2};
            err = var_data{2};
            fitvals = output_fit{2};
            y_vals = fitvals(x_vals);
            xlbl = 'Nanoparticle Concentration';
        elseif strcmp(bound_type,'normBound')
            data = mean_data{1};
            err = var_data{1};
            fitvals = output_fit{1};
            y_vals = fitvals(x2_vals);
            x_vals = x2_vals .* valence;
            lbl = 'Bound TCRs';
            xlbl = ['pMHC Concentration'];
        elseif strcmp(bound_type,'normPhos')
            data = mean_data{2};
            err = var_data{2};
            fitvals = output_fit{2};
            y_vals = fitvals(x2_vals);
            x_vals = x2_vals .* valence;
            lbl = 'Phos TCRs';
            xlbl = ['pMHC Concentration'];
        end
        
        rho30 = x_vals(y_vals > 30);
        
        if isempty(rho30)
            rho30 = 0;
        else
            rho30 = rho30(1);
        end
        
        P30 = [P30, rho30];
    end
    
    x = log10(kd_vals); y = log10(P30)';
    X = [ones(length(x),1),x'];

    b = X \ y;
    
    %cnp = plotCapacity(ax(1), rho_vals, tpc, bound_type, false);
    
    alpha = b(2);
    
    hold(ax,'on');
    scatter(ax, kd_vals, P30, [], colors, 'filled');
    %errorbar(ax, kd_vals, P30, errorMax, 'LineStyle', 'None');
    line = plot(ax, kd_vals, 10^b(1)*kd_vals.^b(2),'Color',colors);
    ylabel(ax,'$P_{15}$');
    xlabel(ax,'$K_D$');
    grid(ax,'on');
    
end