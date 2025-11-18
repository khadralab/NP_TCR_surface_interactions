function [line2] = TCR_NP_Comparison(axes, kd_vals, rho_vals, tpc, valence, colors)
    
    x_vals = linspace(min(rho_vals), max(rho_vals), 1000);
    rho_vals = 10.^rho_vals;
    x_vals = 10.^x_vals;
    x2_vals = x_vals;
    
    koff = kd_vals ./ 10;

    NP50=[];
    NPmax=[];
    TCR50=[];
    TCRmax=[];

    for v = valence
        if tpc == 0
            surface_type = 'uniform';
            file = ['TauLeap/0TPC/v',num2str(v)];
        elseif strcmp(tpc, '10x20')
            surface_type = 'clusters';
            file = ['TauLeap/10x20TPC/v',num2str(v)];
        elseif strcmp(tpc,'5x20')
            surface_type = 'clusters';
            file = ['TauLeap/5x20TPC/v',num2str(v)];
        else
            surface_type = 'clusters';
            file = ['TauLeap/',num2str(tpc),'TPC/v',num2str(v)];
        end
        
        [mean_data, var_data, output_fit] = load_DR(koff, rho_vals, surface_type, file, 1);
        
        bound_type1 = 'bound';
        bound_type2 = 'np';

        % Plot mean Phosphorylated dose-response, Hill fit and shade standard deviation
        if strcmp(bound_type1,'bound')
            data = mean_data{1};
            err = var_data{1};
            fitvals = output_fit{1};
            y_vals = fitvals(x_vals);
            xlbl = 'Nanoparticle Concentration';
        elseif strcmp(bound_type2,'np')
            
        elseif strcmp(bound_type1, 'rescaled')
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
            x_vals = x2_vals .* v;
            lbl = 'Bound TCRs';
            xlbl = ['pMHC Concentration'];
        elseif strcmp(bound_type,'normPhos')
            data = mean_data{2};
            err = var_data{2};
            fitvals = output_fit{2};
            y_vals = fitvals(x2_vals);
            x_vals = x2_vals .* v;
            lbl = 'Phos TCRs';
            xlbl = ['pMHC Concentration'];
        end
        
        TCR50 = [TCR50, fitvals.ec50];
        TCRmax = [TCRmax, fitvals.eMax];
        
        fitvals = output_fit{3};
        
        NP50 = [NP50, fitvals.ec50];
        NPmax = [NPmax, fitvals.eMax];
        
    end
%{
    plot(axes(1), x_vals, valMat')
    xlabel(axes(1), xlbl);
    ylabel(axes(1), '$\Delta$ Bound TCR');
    set(axes(1), 'xscale','log');
    axes(1).ColorOrder = colors;
%}
    
    %p = polyfit(valence,ec50,5);
    
    hold(axes(1),'on');
    scatter(axes(1), TCR50, NP50, [], colors, '*');
    %line1 = plot(axes(1), valence, polyval(p, valence),'Color',colors2);
    set(axes(1), 'yscale','log', 'xscale','log');
    xlabel(axes(1),'TCR EC50');
    ylabel(axes(1),'NP EC50');
    
    %p = polyfit(valence,emax,5);
    
    hold(axes(2), 'on');
    scatter(axes(2), TCRmax, NPmax, [], colors, '*');
    %line2 = plot(axes(2), valence, polyval(p, valence),'Color',colors2);
    xlabel(axes(2), 'TCR EMax');
    ylabel(axes(2), 'NP EMax');
    
end