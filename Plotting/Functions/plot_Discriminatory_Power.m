function y = plot_Discriminatory_Power(ax1, ax2, ax3, threshold_array, valence_vals, surfaces, kp)
    
    rho_vals = linspace(-4,3.5,16);
    rho_vals = 10.^ rho_vals;
    koff_vals = [0.0001, 0.001, 0.01, 0.1];
    %koff_vals = [0.0001, 0.0003, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3];
    np_radius = 20;
    
    %col = [[0,0,0]; [0,0,1]; [0.5,0,0.5]; [1,0,1]; [1,0,0]];
    %colk = [[0,0,0]; [1,0,0]; [0,1,0]];
    %colv = summer(length(valence_vals)+2);

    col = gray(length(valence_vals)+2);
    hol = winter(length(valence_vals)+2);
    sol = autumn(length(valence_vals)+2);
    pol = cool(10);
    mol = summer(length(valence_vals)+2);
    
    linestyle = {'-','--'};
    markerstyle = {'s','d'};
    colorstyle = [col(1,:); hol(1,:); pol(3,:); pol(7,:); sol(1,:)];
        
    %x_vals = linspace(min(log10(rho_vals)), max(log10(rho_vals)), 1000);
    x_vals = linspace(-3.5, max(log10(rho_vals)), 1000);
    x_vals = 10.^x_vals;

    xFit = linspace(min(log10(koff_vals)), max(log10(koff_vals)),1000);
    
    for j = 1:length(surfaces)
        tpc = surfaces(j);

        if ismember(tpc, [3,20])
            %koff_vals = [0.0001, 0.0003, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3];
            koff_vals = [0.0001, 0.001, 0.01, 0.1];
            los=5; loe=8;
            his=1; hie=4;
        elseif ismember(tpc, [1])
            %koff_vals = [0.0001, 0.0003, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3];
            koff_vals = [0.0001, 0.001, 0.01, 0.1];
            los=4; loe=6;
            his=1; hie=3;
        else
            koff_vals = [0.0001, 0.001, 0.01, 0.1];
            los=2; loe=4;
            his=1; hie=3;
        end
        los=2; loe=4;
        his=1; hie=3;
        for v = 1:length(valence_vals)
            valence = valence_vals(v);
            for i = 1:length(koff_vals)
                koff = koff_vals(i);
                [mean_data,var_data,output_fit] = load_DR(koff, kp, rho_vals, tpc, np_radius, valence, 1);
                fitvals = output_fit{2};
                y_vals = fitvals(x_vals);
        
                for k = 1:length(threshold_array)
                    array = 5*x_vals(y_vals > threshold_array(k));
                    if isempty(array)
                        P(i,k) = nan;
                    else
                        P(i,k) = array(1);
                    end

                    if i == length(koff_vals)
                        % Using the first 3 koff for fitting since last one
                        % is low affinity outlier
                        koff_lofit = koff_vals(los:loe);
                        P_lofit = P(los:loe,:);
                        clo(k,:) = polyfit(log10(koff_lofit(~isnan(P_lofit(:,k)))), log10(P_lofit(~isnan(P_lofit(:,k)),k)),1);   
                        ylo(k,:) = polyval(clo(k,:), xFit);

                        koff_hifit = koff_vals(his:hie);
                        P_hifit = P(his:hie,:);
                        chi(k,:) = polyfit(log10(koff_hifit(~isnan(P_hifit(:,k)))), log10(P_hifit(~isnan(P_hifit(:,k)),k)),1);   
                        yhi(k,:) = polyval(chi(k,:), xFit);

                        if ismember(tpc,[1,3,20]) && valence == 5
                            loglog(ax1, 10*koff_vals, P(:,k),'x','Color', colorstyle(j,:),'MarkerSize',10,'HandleVisibility','off')
                            hold(ax1,'on');
                            loglog(ax1, 10*10.^xFit, 10.^ylo(k,:),'LineStyle',':','Color',colorstyle(j,:));
                            loglog(ax1, 10*10.^xFit, 10.^yhi(k,:),"LineStyle","--",'Color',colorstyle(j,:));
                        end         

                        if ~isempty(ismissing(ax2))
                            scatter(ax2, tpc, chi(k,1),72, colorstyle(j,:), 'x');
                            hold(ax2,'on');
                        end

                        if ~isempty(ismissing(ax3))    
                            scatter(ax3, tpc, clo(k,1),72, colorstyle(j,:), 'x');
                            hold(ax3,'on');
                        end
                        
                    end
                end
            end
        end
    end

    %ylim([ax2 ax3], [0. 2.]);
    %ylim(ax4, [-2. 2.]);

    xticks([ax2 ax3],surfaces)
    yticks([ax2 ax3], linspace(0,2,5))
    grid([ax1 ax2 ax3], 'on');

    %ax1.FontSize = 16;
    %ax2.FontSize = 16;
    %ax3.FontSize = 16;
 
    %xlabel(ax1, "$K_D$", 'FontSize',16,'Interpreter','Latex');
    %ylabel(ax1, "pMHC Concentration","FontSize",16,"Interpreter","Latex");
    
    %xlabel([ax2 ax3], "TPC", 'FontSize',16,'Interpreter','Latex')
    %ylabel([ax2 ax3], "Discriminatory Power",'FontSize',16,'Interpreter','latex')

    %xlabel(ax3, "Valence", 'FontSize',16,'Interpreter','Latex')
    %ylabel(ax3, "Discriminatory Power",'FontSize',16,'Interpreter','latex')
    
    title(ax2, 'High Affinity', 'FontSize',16);
    title(ax3, 'Low Affinity', 'FontSize',16);
    
    %l = legend(ax2, string(valence_vals));
    %title(l, "Valences");
end