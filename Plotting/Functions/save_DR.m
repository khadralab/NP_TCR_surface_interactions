function [mean_data, var_data, output_fit] = save_DR(koff, surface_type, valence)
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
    
    bound = []; phos = []; np = [];
        
       
    if strcmp(surface_type,'0')
        for rho = rho_vals
            file = ['LongSims/TauLeap/0TPC/v',num2str(valence),f];
            load([file,'/Uniform_rho',num2str(floor(rho*10000)),'.mat']);
            num_sims = length(homo_bound_tcr);
            data_bt = []; data_pt = []; data_np = [];

            for i=1:num_sims
                bt = homo_bound_tcr{i}(1,end-wind:end);
                pt = homo_phos_tcr{i}(1,end-wind:end);
                bn = homo_bound_np{i}(1,end-wind:end);
                
                bt = bt(~isnan(bt)); pt = pt(~isnan(pt)); bn = bn(~isnan(bn));
                data_bt = [data_bt; bt']; data_pt = [data_pt; pt']; data_np = [data_np; bn'];
            end

            bound = [bound, data_bt]; phos = [phos, data_pt]; np = [np, data_np];

        end
        
    else
        for rho = rho_vals
            file = ['LongSims/TauLeap/',surface_type,'TPC/v',num2str(valence),f];
            load([file,'/Cluster_rho',num2str(floor(rho*10000)),'.mat']);
            num_sims = length(cluster_bound_tcr);
            data_bt = []; data_pt = []; data_np = [];
            
            for i=1:num_sims
                bt = cluster_bound_tcr{i}(1,end-wind:end);
                pt = cluster_phos_tcr{i}(1,end-wind:end);
                bn = cluster_bound_np{i}(1, end-wind:end);
                
                bt = bt(~isnan(bt)); pt = pt(~isnan(pt)); bn = bn(~isnan(bn));
                data_bt = [data_bt; bt']; data_pt = [data_pt; pt']; data_np = [data_np; bn'];
            end
            
            bound = [bound, data_bt]; phos = [phos, data_pt]; np = [np, data_np];

        end
    end
    
    mean_bound = mean(bound,1);
    var_bound = var(bound,0,1);
    mean_phos = mean(phos,1);
    var_phos = var(phos,0,1);
    mean_np = mean(np,1);
    var_np = mean(np,1);
    
    % Fit Hill Function to dose response

    fo = fitoptions('Method','NonlinearLeastSquares',...
                   'Lower',[0,0,0],...
                   'Upper',[Inf,max(rho_vals),10],...
                   'StartPoint',[1 1 1]);
    ft = fittype('eMax*x^n / (x^n + ec50^n)','options',fo);
    
    % Bootstrapping
    B=100;
    for i=1:B
         ind = randi(length(bound), 1, 10);
         
         bootstrap_sample = bound(ind,:);
         
         bootstrap_mean = mean(bootstrap_sample,1);
         
         [fit_bound, gof] = fit(rho_vals', bootstrap_mean' ,ft)
         
         ec50(i) = fit_bound.ec50;
         eMax(i) = fit_bound.eMax;
         n(i) = fit_bound.n;
     end
    
     mean_50 = mean(ec50);
     var_50 = std(ec50);
     
     mean_Max = mean(eMax);
     var_Max = std(eMax);
     
     mean_n = mean(n);
     var_n = std(n);

    [fit_bound,gof] = fit(rho_vals',mean_bound',ft);
    [fit_phos,gof] = fit(rho_vals',mean_phos',ft);
    [fit_np,gof] = fit(rho_vals',mean_np',ft);
    
    out_bound = bound;
    out_phos = phos;
    out_np = np;
    
    mean_data = {out_bound, out_phos, out_np};
    var_data = {var_bound, var_phos, var_np};
    output_fit = {fit_bound, fit_phos, fit_np, var_50, var_Max};
    
    fname = [file,'/dr_curve.mat'];
    save(fname, "mean_data", "var_data", "output_fit")
    
end