function [mean_data, var_data, output_fit] = load_DR(koff, kp, rho_vals, tpc, np_radius, valence, fittrue)
    num_sims = 30;
    wind = 5e1;
    
    x_vals = linspace(min(rho_vals), max(rho_vals), 1000);
    x_vals = 10.^x_vals;
    
    mean_bound = []; var_bound = [];
    mean_phos = []; var_phos = [];
    mean_np = []; var_np = [];

    if tpc == 1
        surface_type = 'uniform';
        file = ['r',num2str(np_radius),'/1TPC/v',num2str(valence)];
    elseif strcmp(tpc, '10x20')
        surface_type = 'clusters';
        file = ['r20/10x20TPC/v',num2str(valence)];
    elseif strcmp(tpc,'5x20')
        surface_type = 'clusters';
        file = ['r20/5x20TPC/v',num2str(valence)];
    elseif strcmp(tpc,'33x3')
        surface_type = 'clusters';
        file = ['r20/33x3TPC/v',num2str(valence)];
    elseif strcmp(tpc,'66x3')
        surface_type = 'clusters';
        file = ['r20/66x3TPC/v',num2str(valence)];
    else
        surface_type = 'clusters';
        file = ['r',num2str(np_radius),'/',num2str(tpc),'TPC/v',num2str(valence)];
    end
    
    
    if mod(koff,0.00003) == 0 && koff < 0.005
        n = abs(log10(koff/3));
        f = ['LongSims/',file,'/koff3e',num2str(n)];        
    elseif koff < 0.005
        n = abs(log10(koff));
        f = ['LongSims/',file,'/koff1e',num2str(n)];
    else
        f = ['LongSims/',file,'/koff',num2str(koff*100)];
    end
    
    bound = []; phos = []; np = [];
        
    if strcmp(surface_type,'clusters')
        for rho = rho_vals
            load([f,'/Cluster_rho',num2str(floor(rho*10000)),'.mat'])
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
        
    else
        for rho = rho_vals
            load([f,'/Uniform_rho',num2str(floor(rho*10000)),'.mat'])
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
    end
    
    % Simple KPR Model
    phos_coef = (kp(1) / (kp(1)+koff))^kp(2);
    phos = phos_coef * bound;

    % Negative Feedback KPR
    %kp_array = kp(1) * (1 - bound / 300);
    %phos_coef = (kp_array ./ (kp_array+koff)).^kp(2);
    %phos = phos_coef .* bound;

    mean_bound = mean(bound,1);
    var_bound = var(bound,0,1);
    mean_phos = mean(phos,1);
    var_phos = var(phos,0,1);
    mean_np = mean(np,1);
    var_np = mean(np,1);
    
    % Fit Hill Function to dose response
    
    if fittrue

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

             [fit_bound, gof] = fit(rho_vals', bootstrap_mean' ,ft);

             ec50(i) = fit_bound.ec50;
             eMax(i) = fit_bound.eMax;
         end

         mean_50 = mean(ec50);
         var_50 = std(ec50);

         mean_Max = mean(eMax);
         var_Max = std(eMax);

        [fit_bound,gof] = fit(rho_vals',mean_bound',ft);
        [fit_phos,gof] = fit(rho_vals',mean_phos',ft);
        [fit_np,gof] = fit(rho_vals',mean_np',ft);

        out_bound = bound;
        out_phos = phos;
        out_np = np;

        mean_data = {out_bound, out_phos, out_np};
        var_data = {var_bound, var_phos, var_np};
        output_fit = {fit_bound, fit_phos, fit_np, var_50, var_Max};
    else
        out_bound = bound;
        out_phos = phos;
        out_np = np;

        mean_data = {out_bound, out_phos, out_np};
        var_data = {var_bound, var_phos, var_np};
        output_fit = 0;
    end
    
end