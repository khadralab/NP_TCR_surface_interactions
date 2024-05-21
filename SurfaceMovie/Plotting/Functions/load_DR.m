function [mean_data, var_data, output_fit] = load_DR(koff, rho_vals, surface_type, file)
    num_sims = 30;
    wind = 5e1;
    
    x_vals = linspace(min(rho_vals), max(rho_vals), 1000);
    x_vals = 10.^x_vals;
    
    mean_bound = []; var_bound = [];
    mean_phos = []; var_phos = [];
    mean_np = []; var_np = [];
    
    if koff < 0.005
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
    
    mean_bound = mean(bound,1);
    var_bound = var(bound,1);
    mean_phos = mean(phos,1);
    var_phos = var(phos,1);
    mean_np = mean(np,1);
    var_np = mean(np,1);
    
    % Fit Hill Function to dose response

    fo = fitoptions('Method','NonlinearLeastSquares',...
                   'Lower',[0,0,0],...
                   'Upper',[Inf,max(rho_vals),10],...
                   'StartPoint',[1 1 1]);
    ft = fittype('eMax*x^n / (x^n + ec50^n)','options',fo);

    [fit_bound,gof] = fit(rho_vals',mean_bound',ft);
    [fit_phos,gof] = fit(rho_vals',mean_phos',ft);
    [fit_np,gof] = fit(rho_vals',mean_np',ft);
    
    out_bound = bound;
    out_phos = phos;
    out_np = np;
    
    mean_data = {out_bound, out_phos, out_np};
    var_data = {var_bound, var_phos, var_np};
    output_fit = {fit_bound, fit_phos, fit_np};
    % Fit Hyperbolic tangent
    %{
    fun = @(par, data) par(1)*tanh(par(2)*data+par(3))+par(4);
    par0 = [50;1;1;50];
    
    options = optimoptions('lsqcurvefit');
    options.MaxFunEvals = 1000;
    options.MaxIterations = 1000;
    
    [par_bound,resnorm,residual,exitflag,output,lambda,jacobian] = ...
        lsqcurvefit(fun,par0,rho_vals',mean_bound');
    
    [par_phos,resnorm,residual,exitflag,output,lambda,jacobian] = ...
        lsqcurvefit(fun,par0,rho_vals',mean_phos');
    %}
end