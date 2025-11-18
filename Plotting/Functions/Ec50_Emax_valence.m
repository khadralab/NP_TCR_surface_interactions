function [n] = Ec50_Emax_valence(ax, kd_vals, rho_vals, surfaces, valence, bound_type, colors)
    x_vals = linspace(min(rho_vals), max(rho_vals), 1000);
    rho_vals = 10.^rho_vals;
    x_vals = 10.^x_vals;
    
    ec50 = []; eMax = []; n=[];
    error50 = []; errorMax = []; errorn = [];

    col = gray(length(valence)+2);
    hol = winter(length(valence)+2);
    sol = autumn(10);
    pol = spring(length(valence)+2);
    mol = summer(10);
    lol = lines(length(kd_vals));
    
    linestyle = {'-','--'};
    markerstyle = {'s','d'};
    %colorstyle = [col(1,:); hol(1,:); mol(1,:); pol(1,:); sol(1,:)];
    
    i=1;
    
    if strcmp(bound_type,'bound')
        k=1;
    elseif strcmp(bound_type,'phos')
        k=2;
    else
        k=3;
    end

    kp = [0.1,2];

    colorstyle = [sol(4,:); sol(6,:); mol(5,:); mol(1,:)];
    markersz = [15,25,30,35,40,50];

    for koff = kd_vals ./ 10
        ec50 = [];
        eMax = [];
        for v = valence
            [mean_data, var_data, output_fit] = load_DR(koff, kp, rho_vals, surfaces, 20, v,1);
            fit_data = output_fit{k};
            ec50 = [ec50,fit_data.ec50];                    % Multiply by v (valence) to normalize to pMHC concentration
            error50 = [error50, output_fit{4}];
            
            
            n = [n,fit_data.n];
            errorn = [errorn, 0];
            
            eMax = [eMax, fit_data.eMax];
            errorMax = [errorMax, output_fit{5}];
        end
        hold(ax,'on');
    
        %pf = polyfit(log10(ec50), eMax, 2);
        %xf = linspace(-3,2,300);
        %yf = polyval(pf, xf);
        
        scatter(ax, ec50, eMax, markersz, colorstyle(i,:), 'filled','LineWidth',1.3,'DisplayName',num2str(koff));
        plot(ax, ec50, eMax, ':', 'Color',colorstyle(i,:),'LineWidth',1.3,'HandleVisibility','off');
        %plot(ax, 10.^xf, yf, 'Color',colors(i,:),'DisplayName',num2str(koff))
        i=i+1;
    end
    
    %l=legend(ax);
    %title(l,"$k_{off}$",'Interpreter','latex');
    %ec50 = log10(ec50);
    
    
end