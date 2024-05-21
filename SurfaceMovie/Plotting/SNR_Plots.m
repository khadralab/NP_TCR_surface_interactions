%% SNR Plots

bound_type = {'bound', 'phos'};
surface_type = {'clusters','uniform'};
rho_vals = linspace(-2,5,15);
rho_vals = 10.^rho_vals;

rho = rho_vals(6);
koff_vals = [0.01, 0.05];

kd_lbls = arrayfun(@num2str, koff_vals ./ 0.1, 'UniformOutput',0);

cols = autumn(length(koff_vals));
hols = winter(length(koff_vals));

f = figure('Color', 'white');
ax(1) = subplot(3,4,1:2); ax(3) = subplot(3,4,3:4); ax(5) = subplot(3,4,9:10);
ax(2) = subplot(3,4,5:6); ax(4) = subplot(3,4,7:8); ax(6) = subplot(3,4,11:12);

for i=5:6
    axis(ax(i),'off');
end

h1 = TCR_Distributions(ax(1), koff_vals, rho, surface_type{1}, bound_type{1},cols);
h2 = TCR_Distributions(ax(2), koff_vals, rho, surface_type{1}, bound_type{2},cols);

h3 = TCR_Distributions(ax(3), koff_vals, rho, surface_type{2}, bound_type{1},hols);
h4 = TCR_Distributions(ax(4), koff_vals, rho, surface_type{2}, bound_type{2},hols);

for i=1:2
    xlabel(ax(2*i-1),'Bound TCRs')
    xlabel(ax(2*i), 'Phos. TCRs')
    ylabel(ax(i),'PDF');
    ylabel(ax(i+2), 'PDF');
    set(ax(i), 'tickdir','out')
    set(ax(i+2),'tickdir','out')
    xlim(ax(i), [0 180]);
    xlim(ax(i+2), [0 60]);
end

% Legend

pos1 = get(ax(5),'Position');
pos2 = get(ax(6),'Position');

l1 = legend(ax(1), h1, kd_lbls,'Position', pos1, 'NumColumns',2);
l2 = legend(ax(3), h3, kd_lbls,'Position',pos2,'NumColumns',2);
title(l1, 'Clusters $K_D$');
title(l2, 'Uniform $K_D$');

set(findall(f,'-property','FontSize'),'FontSize',16)
set(findall(f,'-property','Interpreter'),'Interpreter','Latex')

%%
bound_type = {'bound', 'phos', 'both'};
surface_type = {'clusters','uniform'};
rho_vals = linspace(-2,5,15);
rho_vals = rho_vals(1:10);
rho_vals = 10.^rho_vals;
koff_vals = [0.01, 0.02, 0.05, 0.1];
koff_vals = [0.05];

kd_lbls = arrayfun(@num2str, koff_vals ./ 0.1, 'UniformOutput',0);

cols = autumn(length(koff_vals));
hols = winter(length(koff_vals));

f2 = figure('Color','white');
for i=1:1
    ax(i) = subplot(1,1,i);
    set(ax(i),'xscale','log','yscale','log','tickdir','out')
    xlabel(ax(i),'NP Concentration')
    grid(ax(i),'on');
end

h1 = Plotting_SNR(ax(1), koff_vals, rho_vals, surface_type{1}, bound_type{1}, hols);
h2 = Plotting_SNR(ax(1), koff_vals, rho_vals, surface_type{2}, bound_type{1}, cols);

h1.MarkerSize=15;
h2.MarkerSize=15;

%Plotting_SNR(ax(1), koff_vals, rho_vals, surface_type{1}, bound_type{2}, cols);
%Plotting_SNR(ax(2), koff_vals, rho_vals, surface_type{2}, bound_type{2}, cols);
%Plotting_SNR_Ratio(ax(3), koff_vals, rho_vals, surface_type{1}, cols)
%Plotting_SNR_Ratio(ax(4), koff_vals, rho_vals, surface_type{2}, hols)

for i=1:1
    ylabel(ax(i), 'SNR')
    ylim(ax(i), [0 1e3]);

    %ylabel(ax(i+2), 'Ratio SNR')
    %ylabel(ax(i+4), 'Phos / Bound')
end

%ax(7) = subplot(427); ax(8) = subplot(428);
%axis(ax(7),'off'); axis(ax(8),'off');

%pos1 = get(ax(7),'Position');
%pos2 = get(ax(8),'Position');

legend(ax(1), [h1 h2], {'Clusters', 'Uniform'});

%l1 = legend(ax(1), h1, bound_type,'Position',pos1, 'NumColumns',2);
%l2 = legend(ax(2), h2, kd_lbls,'Position',pos2,'NumColumns',2);
%title(ax(1), 'Clusters');
%title(ax(2), 'Uniform');

set(findall(f2,'-property','FontSize'),'FontSize',24)
set(findall(f2,'-property','Interpreter'),'Interpreter','Latex')

%% Local Functions
%% Plotting Histograms

function [H, average, variance] = TCR_Distributions(ax, koff_vals, rho, surface_type, bound_type, cols)
    hold(ax,'on');
    
    i=1;
    for koff = koff_vals
        [tcr_dist] = load_DR(koff, rho, surface_type, bound_type);
        
        start = floor(0.75*length(tcr_dist));
        
        average(i) = mean(tcr_dist(:,start:end),'all');
        variance(i) = var(tcr_dist(:,start:end),0,'all');
    
        H(i) = histogram(ax, tcr_dist(:,start:end),'FaceAlpha',0.2,'EdgeAlpha',0.2,'FaceColor',cols(i,:),'Normalization','pdf','BinMethod','integers');
        xline(ax,average(i),'--k',{[num2str(floor(average(i)))]},'Linewidth',1.5, 'LabelOrientation', 'Horizontal');
        i=i+1;
    end
    
end

%% Plotting SNR Profiles

function [y] = Plotting_SNR(ax, koff_vals, rho_vals, surface_type, bound_type, cols)

    hold(ax, 'on')
    i=1;
    for koff = koff_vals
        x=[];
        for rho = rho_vals
            [snr, mean_snr, var_snr] = Bootstrap_SNR(koff, rho, surface_type, bound_type);
            x = [x, mean_snr];
        end
        y(i) = plot(ax, rho_vals, x,'-*','Color',cols(i,:));
        
        i=i+1;
    end
    
end

function [line] = Plotting_SNR_Ratio(ax, koff_vals, rho_vals, surface_type, cols)

    hold(ax, 'on')
    i=1;
    for koff = koff_vals
        
        x=[]; y =[];
        
        for rho = rho_vals
            [snr, mean_snr, var_snr] = Bootstrap_SNR(koff, rho, surface_type, 'phos');
            x = [x, mean_snr];
            [snr, mean_snr, var_snr] = Bootstrap_SNR(koff, rho, surface_type, 'bound');
            y = [y, mean_snr];
        end
        line(i) = plot(ax, rho_vals, x./y,'-*','Color',cols(i,:));
        
        i=i+1;
    end
    
end
%% Loading Data

function [tcr_dist] = load_DR(koff, rho, surface_type, bound_type)
        
    if strcmp(surface_type,'clusters')
        load(['LongSims/TauLeap/koff',num2str(koff*100),'/Cluster_rho',num2str(floor(rho*100)),'.mat']);
        bound_tcr = cluster_bound_tcr;
        phos_tcr = cluster_phos_tcr;

    else
        load(['LongSims/TauLeap/koff',num2str(koff*100),'/Uniform_rho',num2str(floor(rho*100)),'.mat']);
        bound_tcr = homo_bound_tcr;
        phos_tcr = homo_phos_tcr;
    end
    
    time_points = 1e5;
    end_time = 1e6;
    time_array = linspace(0,end_time,time_points);
    
    for i = 1:size(bound_tcr,2)
        y1(i,:) = interp1(bound_tcr{i}(2,:), bound_tcr{i}(1,:), time_array);
        y2(i,:) = interp1(phos_tcr{i}(2,:), phos_tcr{i}(1,:), time_array);
        
    end
    
    if strcmp(bound_type,'phos')
        tcr_dist = y2;
    else
        tcr_dist = y1;
    end
    
end

%% Bootstrap Method

function [snr, mean_snr, var_snr] = Bootstrap_SNR(koff, rho, surface_type, bound_type)
    tcr_dist = load_DR(koff, rho, surface_type, bound_type);
    
    n = size(tcr_dist,1);
    m = size(tcr_dist,2);
    k = floor(0.2*n);
    start = floor(0.5*m);  % Select the last 25% of time series
    
    snr = zeros(1,n);
    
    i=1;
    
    while i <= n
        ind = randsample(n, k);
        x = tcr_dist(ind, start:end);
        
        snr(i) = mean(x.^2,'all') ./ var(x,0,'all','omitnan');
        
        i=i+1;
    end
    
    mean_snr = mean(snr); var_snr = var(snr);
end