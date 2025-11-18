function [H, average, variance] = TCR_Distributions(ax, koff_vals, rho_vals, valence, tpc, bound_type, cols)
    hold(ax,'on');
    
    for i=1:length(koff_vals)
        koff = koff_vals(i); rho = 10^rho_vals(i); v = valence(i);
        
        tcr_dist = load_data(koff, rho, v, tpc, 'tcr');
        
        start = floor(0.75*length(tcr_dist));
        
        average(i) = mean(tcr_dist(:,start:end),'all');
        variance(i) = var(tcr_dist(:,start:end),0,'all');
    
        H(i) = histogram(ax, tcr_dist(:,start:end),'FaceAlpha',0.8,'EdgeAlpha',0.2,'FaceColor',cols(i,:),'Normalization','pdf','BinMethod','integers');
        %xline(ax,average(i),'--k','','Linewidth',1.5, 'LabelOrientation', 'Horizontal');
        i=i+1;
    end
    
    
    if strcmp(bound_type,'phos')
        xlbl = 'Phos. TCRs';
    elseif strcmp(bound_type,'np')
        xlbl = 'Bound NPs';
    else
        xlbl = 'Bound TCRs';
    end
    xlabel(ax, xlbl);
    
end

%% Loading Data

function [tcr_dist] = load_data(koff, rho, v, tpc, bound_type)

    file = ['LongSims/r20/',num2str(tpc),'TPC/v',num2str(v),'/'];
    
    if tpc == 1
        file = ['LongSims/r20/1TPC/v',num2str(v),'/'];
    end
        
    if mod(koff,0.00003) == 0 && koff < 0.005
        n = abs(log10(koff/3));
        f = [file,'koff3e',num2str(n)];        
    elseif koff < 0.005
        n = abs(log10(koff));
        f = [file,'koff1e',num2str(n)];
    else
        f = [file,'koff',num2str(koff*100)];
    end
        
    if tpc == 1
        load([f,'/Uniform_rho',num2str(floor(rho*10000)),'.mat']);
        bound_tcr = homo_bound_tcr;
        bound_np = homo_bound_np;
        phos_tcr = homo_phos_tcr;

    else
        load([f,'/Cluster_rho',num2str(floor(rho*10000)),'.mat']);
        bound_tcr = cluster_bound_tcr;
        bound_np = cluster_bound_np;
        phos_tcr = cluster_phos_tcr;
    end
    
    time_points = 1e5;
    end_time = 1e6;
    time_array = linspace(0,end_time,time_points);
    
    for i = 1:size(bound_tcr,2)
        y1(i,:) = interp1(bound_tcr{i}(2,:), bound_tcr{i}(1,:), time_array);
        y2(i,:) = interp1(phos_tcr{i}(2,:), phos_tcr{i}(1,:), time_array);
        y3(i,:) = interp1(bound_np{i}(2,:), bound_np{i}(1,:), time_array);
        
    end
    
    if strcmp(bound_type,'phos')
        tcr_dist = y2;
    elseif strcmp(bound_type,'np')
        tcr_dist = y3;
    else
        tcr_dist = y1;
    end
    
end
