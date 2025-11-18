function [mean_data, var_data, output_fit] = load_DR2(koff, surface_type, valence)
    num_sims = 30;
    wind = 5e1;
    
    rho_vals = linspace(-4,3.5,16);
    rho_vals = 10.^rho_vals

    x_vals = linspace(min(rho_vals), max(rho_vals), 1000);
    x_vals = 10.^x_vals;
    
    
    mean_bound = []; var_bound = [];
    mean_phos = []; var_phos = [];
    mean_np = []; var_np = [];
    
    if mod(koff,0.00003) == 0 && koff < 0.005
        n = abs(log10(koff/3));
        f = ['/koff3e',num2str(n)];        
    elseif koff < 0.005
        n = abs(log10(koff));
        f = ['/koff1e',num2str(n)];
    else
        f = ['/koff',num2str(koff*100)];
    end
    
    if strcmp(surface_type, '0')
        file = ['LongSims/TauLeap/',surface_type,'TPC/v',num2str(valence),f];
        
    else
        file = ['LongSims/TauLeap/',surface_type,'TPC/v',num2str(valence),f];
    end
    
    load([file,'/dr_curve.mat']);