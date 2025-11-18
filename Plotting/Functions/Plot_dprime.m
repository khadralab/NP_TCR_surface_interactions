function [dp, x] = Plot_dprime(ax, koff_vals, rho_vals, valence, tpc, bound_type, cols, arg)
    if arg
        x(:,1) = rho_vals(1:end-1)';
        x(:,2) = rho_vals(2:end)';
        
        xlbl = 'Nanoparticle Concentration';
        
        dp = zeros(length(x),1);
        
        for i=1:length(x)
            dp(i) = dprime(koff_vals, x(i,:), valence, tpc, bound_type);
        end
    else
        x(:,1) = koff_vals(1:end-1)';
        x(:,2) = koff_vals(2:end)';
        
        xlbl = 'Koff';
        
        dp = zeros(length(x),1);
        
        for i=1:length(x)
            dp(i) = dprime(x(i,:), rho_vals, valence, tpc, bound_type);
        end
    end

    %plot(ax,x(:,1),dp,':x','Color',cols,'Linewidth',1.5)
    colororder(cols);
    h = plot(ax,rho_vals(1:end-1)',dp,':x','Linewidth',1.5);
    ylabel(ax,' d-prime')
    xlabel(ax, xlbl);
    
end

function [dp] = dprime(koff_vals, rho_vals, valence, tpc, bound_type)

    mu = [0,0];
    variance = [0,0];
    
    for i=1:2
       tcr_dist = load_DR(koff_vals(i), rho_vals(i), valence(i), tpc, bound_type);
       start = floor(0.99*length(tcr_dist));
       mu(i) = mean(tcr_dist(:,start:end),'all');
       variance(i) = var(tcr_dist(:, start:end),0,'all');
    end
    
    dp = abs(mu(1) - mu(2)) / sqrt(0.5*(variance(1) + variance(2)));
end

function [tcr_dist] = load_DR(koff, rho, v, tpc, bound_type)

    file = ['LongSims/TauLeap/',num2str(tpc),'TPC/v',num2str(v),'/'];
    
    if tpc == 0 
        file = ['LongSims/TauLeap/20TPC/v',num2str(v),'/'];
    elseif strcmp(tpc,'20x10')
        file = ['LongSims/TauLeap/20x10TPC/v',num2str(v),'/'];
    elseif strcmp(tpc,'20x5')
        file = ['LongSims/TauLeap/20x5TPC/v',num2str(v),'/'];
    end
        
    if mod(koff,0.00003) == 0 && koff < 0.005
        n = abs(log10(koff/3));
        f = [file,'koff3e',num2str(n)];        
    elseif koff < 0.005
        n = abs(log10(koff));
        f = [file,'/koff1e',num2str(n)];
    else
        f = [file,'/koff',num2str(koff*100)];
    end
        
    if tpc == 0
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