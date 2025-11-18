function [bound_tcr, bound_np] = load_BoundTCR_dist(koff, np_rho, np_radius, np_valence, TPC)

    bound_tcr = []; phos = []; bound_np = [];
    % Retrieve file name should be a seperate function

    if TPC == 1
        surface_type = 'uniform';
    else
        surface_type = 'clusters';
    end

    file = ['LongSims/r',num2str(np_radius),'/',num2str(TPC),'TPC/v',num2str(np_valence)];

    if mod(koff,0.00003) == 0 && koff < 0.005
        n = abs(log10(koff/3));
        f = [file,'/koff3e',num2str(n)];        
    elseif koff < 0.005
        n = abs(log10(koff));
        f = [file,'/koff1e',num2str(n)];
    else
        f = [file,'/koff',num2str(koff*100)];
    end
       
    wind = 50;
    if strcmp(surface_type,'clusters')
        load([f,'/Cluster_rho',num2str(floor(np_rho*10000)),'.mat'])
        num_sims = length(cluster_bound_tcr);
        data_bt = []; data_pt = []; data_np = [];
        
        for i=1:num_sims
            bt = cluster_bound_tcr{i}(1,end-wind:end);
            pt = cluster_phos_tcr{i}(1,end-wind:end);
            bn = cluster_bound_np{i}(1, end-wind:end);
            bt = bt(~isnan(bt)); pt = pt(~isnan(pt)); bn = bn(~isnan(bn));
            data_bt = [data_bt; bt']; data_pt = [data_pt; pt']; data_np = [data_np; bn'];
        end
        
        bound_tcr = data_bt; phos = data_pt; bound_np = data_np;
        
    else
        load([f,'/Uniform_rho',num2str(floor(np_rho*10000)),'.mat'])
        num_sims = length(homo_bound_tcr);
        data_bt = []; data_pt = []; data_np = [];

        for i=1:num_sims
            bt = homo_bound_tcr{i}(1,end-wind:end);
            pt = homo_phos_tcr{i}(1,end-wind:end);
            bn = homo_bound_np{i}(1,end-wind:end);
            
            bt = bt(~isnan(bt)); pt = pt(~isnan(pt)); bn = bn(~isnan(bn));
            data_bt = [data_bt; bt']; data_pt = [data_pt; pt']; data_np = [data_np; bn'];
        end

        bound_tcr = data_bt; phos = data_pt; bound_np = data_np;
    end
end