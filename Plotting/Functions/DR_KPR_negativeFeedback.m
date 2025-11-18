function [scat] = DR_KPR_negativeFeedback(ax, kd_vals, kp, rho_vals, tpc, valence, colors, xlb, ylb, np_radius, feedbacktrue)
    x_vals = linspace(min(rho_vals), max(rho_vals), 1000);
    rho_vals = 10.^rho_vals;
    x_vals = 10.^x_vals;
    fun = @(par, data) par(1)*tanh(par(2)*data+par(3))+par(4);
    xlbl = ['NP Concentration'];
    
    i=1;

    kp=[kp, 0.07, 5];
    ST=0.14;
    Cs=200;
    
    x2_vals = valence*x_vals;
    np_rho = rho_vals;

    for koff = kd_vals ./ 10
        [mean_data, var_data, output_fit] = load_DR(koff, kp, np_rho, tpc, np_radius, valence, 1);
       
        % Plot mean Phosphorylated dose-response, Hill fit and shade standard deviation
        data = mean_data{1};
        err = var_data{1};
        fitvals = output_fit{1};
        y_vals = fitvals(x_vals);
        for k=1:length(y_vals)
            y_vals(k) = negative_feedback_KPR(y_vals(k), koff, kp, ST, Cs, feedbacktrue);
        end

        lbl = 'Phos TCRs';
                
        hold(ax, 'on');
        %axes(ax);
        %boxplot(ax,data,rho_vals,'MedianStyle','target','Positions',rho_vals,'Labels',rho_vals,'colors', colors(i,:),'BoxStyle','filled','Symbol','','Widths', 0.5 * rho_vals);
        scat(i) = plot(ax,x2_vals, y_vals, 'Color', colors(i,:), 'Linewidth', 1.3);
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