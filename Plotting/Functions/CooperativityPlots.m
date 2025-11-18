function [alpha, line] = CooperativityPlots(ax, kd_vals, kp, rho_vals, tpc, valence, np_radius, bound_type, colors, xlb, ylb, num_plts)
    x_vals = linspace(min(rho_vals), max(rho_vals), 1000);
    rho_vals = 10.^rho_vals;
    x_vals = 10.^x_vals;
    x2_vals = x_vals;
    
    %{
    if tpc == 1
        surface_type = 'uniform';
        file = [fname,'/0TPC/v',num2str(valence)];
    else
        surface_type = 'clusters';
        file = [fname,'/',num2str(tpc),'TPC/v',num2str(valence)];
    end
    %}

    ec50 = [];
    error50 = [];
    eMax = []; y_vals = [];
    errorMax = [];
    
    i=1;
    
    for koff = kd_vals ./ 10
        [mean_data, var_data, output_fit] = load_DR(koff, kp, rho_vals, tpc, np_radius, valence, 1);
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
        ec50 = [ec50,fitvals.ec50];
        error50 = [error50, output_fit{4}];
        
        eMax = [eMax,fitvals.eMax];
        errorMax = [errorMax, output_fit{5}];
    end
    %{
    if strcmp(bound_type, 'xPMHC')
        ec50 = ec50 .* valence;
    end
    %}
    
    alpha = [];
    
    if num_plts == 1
        x = log10(kd_vals); y = log10(ec50)';
        X = [ones(length(x),1),x'];

        b = X \ y;
        
        p = polyfit(log10(kd_vals),log10(ec50),1);

        alpha = b(2);

        hold(ax,'on');
        %scatter(ax, kd_vals, ec50, [], colors, 'filled');
        errorbar(ax, kd_vals, ec50, error50, '-s', 'Color', colors, 'LineStyle', 'None');
        line = plot(ax, kd_vals, 10.^polyval(p, log10(kd_vals)),'Color',colors);
        %title(ax, ['Slope = ',num2str(alpha)]);
        ylbl = 'EC50 (log$_{10}$)';
        xlbl = '$K_D$ (log$_{10}$ a.u.)';
        grid(ax,'on');
        
    elseif num_plts == 2
        x = log10(kd_vals); y = log10(eMax)';
        X = [ones(length(x),1),x'];
        
        hill_fit = @(b,x) b(1).*b(2) ./ (b(2) + x);
        tanh_fit = @(b,x) b(1).*tanh(b(2).*x + b(3)) + b(4);
        
        %b0 = [-200 1 5 200];
        b0 = [-100 0.1 3.3 100];       % Original
        if valence == 1
            b0 = [-50 0.01 0.01 50];
        end
        lb = [-1000 0 0 -1000];
        ub = [0 2 5 1000];
        
        options = optimoptions('lsqcurvefit','MaxFunctionEvaluations',5000, 'MaxIterations',1000,'OptimalityTolerance',1e-8,'FunctionTolerance',1e-8);

        [B, resnorm, residual, exitflag, output] = lsqcurvefit(tanh_fit, b0, log10(kd_vals),eMax, [], [], options);
        xx = linspace(log10(min(kd_vals)), log10(max(kd_vals)), 1000);

        hold(ax,'on');
        errorbar(ax, kd_vals, eMax, errorMax, '-s', 'Color', colors, 'LineStyle', 'None');
        line = plot(ax, 10.^xx, tanh_fit(B, xx),'Color',colors);
        ylbl = 'EMax';
        xlbl = '$K_D$ (log$_{10}$ a.u.)';
        grid(ax,'on');
    
    elseif num_plts == 3
        p = polyfit(eMax, ec50, 5);
        
        hold(ax,'on');
        errorbar(ax, eMax, ec50, error50, errorMax, '-s', 'Color', colors, 'LineStyle', 'None');
        line = plot(ax, eMax, polyval(p, eMax),'Color',colors);
        
        ylbl = 'EC50';
        xlbl = 'EMax';
        grid(ax,'on');
    end
    
    if xlb
        xlabel(ax, xlbl, 'FontSize', 24)
    end
    
    if ylb
        ylabel(ax,ylbl, 'FontSize', 24)
    end
    
end