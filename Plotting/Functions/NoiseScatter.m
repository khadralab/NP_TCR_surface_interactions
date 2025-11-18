function [y] = NoiseScatter(ax, kd_vals, rho_vals, tpc, valence, bound_type, colors)

    x_vals = linspace(min(rho_vals), max(rho_vals), 1000);
    rho_vals = 10.^rho_vals;
    x_vals = 10.^x_vals;

    dr_vals = [];

    if tpc == 0
        surface_type = 'uniform';
        file = ['TauLeap/0TPC/v',num2str(valence)];
    elseif strcmp(tpc, '10x20')
        surface_type = 'clusters';
        file = ['TauLeap/10x20TPC/v',num2str(valence)];
    elseif strcmp(tpc,'5x20')
        surface_type = 'clusters';
        file = ['TauLeap/5x20TPC/v',num2str(valence)];
    else
        surface_type = 'clusters';
        file = ['TauLeap/',num2str(tpc),'TPC/v',num2str(valence)];
    end

    i=1;

    x2_vals = x_vals;
    np_rho = rho_vals;

    hold(ax, 'on');
    grid(ax, 'on');
    set(ax, 'box','off');

    for koff = kd_vals ./ 10
        [mean_data, var_data, output_fit] = load_DR(koff, np_rho, surface_type, file, 1);

        % Plot mean Phosphorylated dose-response, Hill fit and shade standard deviation
        if strcmp(bound_type,'bound')
            data = mean_data{1};
            err = var_data{1};
            fitvals = output_fit{1};
            y_vals = fitvals(x_vals);
            lbl = 'Bound TCRs';
        elseif strcmp(bound_type, 'rescaled')
            fitvals = output_fit{1};
            data = mean_data{1} ./ fitvals.eMax;
            y_vals = fitvals(x2_vals) ./ fitvals.eMax;
            x_vals = x2_vals ./ fitvals.ec50;
            rho_vals = np_rho ./ fitvals.ec50;
            lbl = 'Relative Activation';
        elseif strcmp(bound_type,'phos')
            data = mean_data{2};
            err = var_data{2};
            fitvals = output_fit{2};
            y_vals = fitvals(x_vals);
            lbl = 'Phosphorylated TCRs';
        elseif strcmp(bound_type,'np')
            data = mean_data{3};
            err = var_data{3};
            fitvals = output_fit{3};
            y_vals = fitvals(x_vals);
            lbl = 'Bound NPs';
        elseif strcmp(bound_type,'covered')
            data = mean_data{1} ./ mean_data{3};
            y_vals = output_fit{1}(x_vals) ./ output_fit{3}(x_vals);
            lbl = 'Bound TCR per NP';             
        elseif strcmp(bound_type,'fracphos')
            data = mean_data{2} ./ mean_data{1};
            y_vals = output_fit{2}(x_vals) ./ output_fit{1}(x_vals);
            lbl = 'Fraction of Phosphorylated TCRs';
        elseif strcmp(bound_type,'xPMHC')
            data = mean_data{1};
            err = var_data{1};
            fitvals = output_fit{1};
            y_vals = fitvals(x2_vals);
            x_vals = x2_vals .* valence;
            rho_vals = np_rho .* valence;
            lbl = 'Bound TCRs';
            xlbl = ['pMHC Concentration'];
        end

        scatter(ax, mean(data), std(data), 25, 'filled')
    end
    
    ax.ColorOrder = colors;
    xlabel(ax, 'Mean')
    ylabel(ax, 'Variance')
    title(ax, [num2str(tpc), 'TPC, valence = ', num2str(valence),':',bound_type])
end
