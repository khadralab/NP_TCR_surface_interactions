function [scat] = DoseResponsePlots2(ax, rho_vals, yvals, yerr, output_fit, colors, xlb, ylb)
        
        x_vals = linspace(min(rho_vals), max(rho_vals), 1000);
        rho_vals = 10.^rho_vals;
        x_vals = 10.^x_vals;

        hold(ax, 'on');
        axes(ax);
        boxplot(yvals,rho_vals,'MedianStyle','target','Positions',rho_vals,'Labels',rho_vals,'colors', colors(i,:),'BoxStyle','filled','Symbol','','Widths', 0.5 * rho_vals);
        scat(i) = plot(ax,x_vals, output_fit{1}(xvals), 'Color', colors(i,:), 'Linewidth', 1.3);
        %shade(ax, rho_vals, data_bound+std_bound, 'w', rho_vals, max(0,data_bound-std_bound),'w', 'FillType',[1 2;2 1],'FillColor', colors(i,:),'FillAlpha',0.1, 'HandleVisibility', 'off');
        grid(ax, 'on');
        
        if ylb
            ylabel(ax, lbl);
        end
        
        if xlb
            xlabel(ax, xlbl);
        end

end