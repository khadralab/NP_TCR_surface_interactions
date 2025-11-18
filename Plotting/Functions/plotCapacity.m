function [carrying_capacity] = plotCapacity(ax, tcrs_per_cluster, rNP, bound_type, cols, plt)
    rho = 10^3.5;
    wind = 1e3;
    carrying_capacity = [];
    %addpath('/home/louis/Desktop/NP-Revival/SurfaceMovie/LongSims');
    
    for tpc = tcrs_per_cluster
        data_bt = []; data_pt = []; data_np = [];
            
        if tpc == 1
            surface_type = 'uniform';
            file = ['LongSims/r',num2str(rNP),'/1TPC/v5/koff0'];
            load([file,'/Uniform_rho',num2str(floor(rho*10000)),'.mat'])
            num_sims = length(homo_bound_tcr);
            
            for i=1:num_sims
                bt = homo_bound_tcr{i}(1,end-wind:end);
                pt = homo_phos_tcr{i}(1,end-wind:end);
                bn = homo_bound_np{i}(1,end-wind:end);
                
                bt = bt(~isnan(bt)); pt = pt(~isnan(pt)); bn = bn(~isnan(bn));
                data_bt = [data_bt; bt']; data_pt = [data_pt; pt']; data_np = [data_np; bn'];
            end
        else
            surface_type = 'clusters';
            file = ['LongSims/r',num2str(rNP),'/',num2str(tpc),'TPC/v5/koff0'];
            
            load([file,'/Cluster_rho',num2str(floor(rho*10000)),'.mat'])
            num_sims = length(cluster_bound_tcr);
            
            for i=1:num_sims
                bt = cluster_bound_tcr{i}(1,end-wind:end);
                pt = cluster_phos_tcr{i}(1,end-wind:end);
                bn = cluster_bound_np{i}(1, end-wind:end);
                
                bt = bt(~isnan(bt)); pt = pt(~isnan(pt)); bn = bn(~isnan(bn));
                data_bt = [data_bt; bt']; data_pt = [data_pt; pt']; data_np = [data_np; bn'];
            end
        end

        if strcmp(bound_type,'bound')
            data = mean(data_bt,'all');
            lbl = 'Total covered TCRs';
        elseif strcmp(bound_type,'np')
            data = mean(data_np,'all');
            lbl = 'NP surface capacity.';
        elseif strcmp(bound_type,'covered')
            data = mean(data_bt ./ data_np, 'all');
            lbl = 'Covered TCRs per NP';   
        elseif strcmp(bound_type,'nc')
            data = mean(data_np ./ (300 / tpc),'all');
            lbl = 'Average bound NPs per nanocluster';
        end
        
        carrying_capacity = [carrying_capacity, data];
        
    end
    
    if plt == true
        plot(ax,tcrs_per_cluster, carrying_capacity,':sq','Color',cols,'DisplayName',[num2str(rNP)],'LineWidth',2);
        xlabel(ax, 'TCRs per Cluster')
        ylabel(ax, lbl)
    end
end