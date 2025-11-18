function [scat] = DoseResponsePlots(ax, kd_vals, kp, rho_vals, tpc, valence, bound_type, colors, xlb, ylb, np_radius)
    x_vals = linspace(min(rho_vals), max(rho_vals), 1000);
    rho_vals = 10.^rho_vals;
    x_vals = 10.^x_vals;
    fun = @(par, data) par(1)*tanh(par(2)*data+par(3))+par(4);
    xlbl = ['NP Concentration'];
    
    i=1;
    
    x2_vals = x_vals;
    np_rho = rho_vals;

    for koff = kd_vals ./ 10
        [mean_data, var_data, output_fit] = load_DR(koff, kp, np_rho, tpc, np_radius, valence, 1);
       
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
        elseif strcmp(bound_type,'normBound')
            data = mean_data{1};
            err = var_data{1};
            fitvals = output_fit{1};
            y_vals = fitvals(x2_vals);
            x_vals = x2_vals .* valence;
            rho_vals = np_rho .* valence;
            lbl = 'Bound TCRs';
            xlbl = ['pMHC Concentration'];
        elseif strcmp(bound_type,'normPhos')
            data = mean_data{2};
            err = var_data{2};
            fitvals = output_fit{2};
            y_vals = fitvals(x2_vals);
            x_vals = x2_vals .* valence;
            rho_vals = np_rho .* valence;
            lbl = 'Phos TCRs';
            xlbl = ['pMHC Concentration'];
        end
                
        hold(ax, 'on');
        %axes(ax);
        %boxplot(ax,data,rho_vals,'MedianStyle','target','Positions',rho_vals,'Labels',rho_vals,'colors', colors(i,:),'BoxStyle','filled','Symbol','','Widths', 0.5 * rho_vals);
        scat(i) = plot(ax,x_vals, y_vals, 'Color', colors(i,:), 'Linewidth', 1.3);
        %shade(ax, rho_vals, data_bound+std_bound, 'w', rho_vals, max(0,data_bound-std_bound),'w', 'FillType',[1 2;2 1],'FillColor', colors(i,:),'FillAlpha',0.1, 'HandleVisibility', 'off');
        grid(ax, 'on');
        
        if ylb
            ylabel(ax, lbl);
        end
        
        if xlb
            xlabel(ax, xlbl);
        end
        
        i=i+1;
    end
    
end