function [scat] = DoseResponsePlots(ax, kd_vals, rho_vals, surface_type, bound_type, colors, file)
    x_vals = linspace(min(rho_vals), max(rho_vals), 1000);
    rho_vals = 10.^rho_vals;
    x_vals = 10.^x_vals;
    fun = @(par, data) par(1)*tanh(par(2)*data+par(3))+par(4);
    
    i=1;

    for koff = kd_vals ./ 10
        [mean_data, var_data, output_fit] = load_DR(koff, rho_vals, surface_type, file);
       
        % Plot mean Phosphorylated dose-response, Hill fit and shade standard deviation
        if strcmp(bound_type,'bound')
            data = mean_data{1};
            err = var_data{1};
            fitvals = output_fit{1};
            y_vals = fitvals(x_vals);
            lbl = 'Bound TCRs';
        elseif strcmp(bound_type,'phos')
            data = mean_data{2};
            err = var_data{2};
            fitvals = output_fit{2};
            y_vals = fitvals(x-Vals);
            lbl = 'Phos TCRs';
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
        end
        
        
        hold(ax, 'on');
        axes(ax);
        boxplot(data,rho_vals,'MedianStyle','target','Positions',rho_vals,'Labels',rho_vals,'colors', colors(i,:),'BoxStyle','filled','Symbol','','Widths', 0.5 * rho_vals);
        scat(i) = plot(ax,x_vals, y_vals, 'Color', colors(i,:));
        %shade(ax, rho_vals, data_bound+std_bound, 'w', rho_vals, max(0,data_bound-std_bound),'w', 'FillType',[1 2;2 1],'FillColor', colors(i,:),'FillAlpha',0.1, 'HandleVisibility', 'off');
        grid(ax, 'on');
        ylabel(ax, lbl);
        xlabel(ax,'NP Concentration');
        
        i=i+1;
    end
    
end