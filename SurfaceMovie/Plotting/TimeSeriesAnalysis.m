%% Time Series Analysis
clear all

file = ['LongSims/TauLeap/Cluster100/koff10'];
rho_vals = linspace(-4,3.5,16);
rho_vals = 10.^rho_vals;

rho = rho_vals(16);

load([file,'/Cluster_rho',num2str(floor(10000*rho)),'.mat']);

%load([file,'/Uniform_rho',num2str(floor(100*rho)),'.mat']);
%% Create Figures for Time Series

f1 = figure('Color','white');
ax(1) = subplot(231); ax(2) = subplot(232); ax(3) = subplot(233);
ax(4) = subplot(234); ax(5) = subplot(235); ax(6) = subplot(236);

for i=1:3 
    % Axis scale for Time Series
    set(ax(i),'tickdir','out')
    
    % Axis scales for Difference Series
    set(ax(i+3),'tickdir','out')
    %
end


TimeSeries(ax(1), cluster_bound_np);
TimeSeries(ax(2), cluster_bound_tcr);
TimeSeries(ax(3), cluster_phos_tcr);

ylabel(ax(1), 'Bound NPs')
ylabel(ax(2), 'Bound TCRs')
ylabel(ax(3), 'Phos TCRs')

title(ax(1:3), 'Cluster Surface');
%{
TimeSeries(ax(4), homo_bound_np);
TimeSeries(ax(5), homo_bound_tcr);
TimeSeries(ax(6), homo_phos_tcr);

ylabel(ax(4), 'Bound NPs')
ylabel(ax(5), 'Bound TCRs')
ylabel(ax(6), 'Phos TCRs')

title(ax(4:6), 'Uniform Surface')
%}
%% Sliding Window

f2 = figure('Color','white');
ax(1) = subplot(231); ax(2) = subplot(232); ax(3) = subplot(233);
ax(4) = subplot(234); ax(5) = subplot(235); ax(6) = subplot(236);

for i=1:3 
    % Axis scale for Cluster Series
    set(ax(i),'tickdir','out')
    
    % Axis scales for Uniform Series
    set(ax(i+3),'tickdir','out')
    ylim(ax([i i+3]), [0 120])

end

rho_vals = [1,100,31622];
for i = 1:3
    rho= rho_vals(i);
    
    load([file,'/Cluster_rho',num2str(rho),'.mat']);
    load([file,'/Uniform_rho',num2str(rho),'.mat']);
    
    slidingWindow(ax(i), cluster_phos_tcr, 1e3);
    title(ax(i), ['$\rho$ = ', num2str(rho)], 'Interpreter','Latex')
    ylabel(ax(1), 'Cluster Phos TCR')

    slidingWindow(ax(i+3), homo_phos_tcr, 1e3);
    ylabel(ax(4), 'Uniform Phos TCR')
end

%% Convergence of Low Concentrations
rho=1000;
load([file,'/Cluster_rho',num2str(rho),'.mat']);
load([file,'/Uniform_rho',num2str(rho),'.mat']);

f3 = figure('Color','white');
ax(1) = subplot(331); ax(2) = subplot(332); ax(3) = subplot(333);
ax(4) = subplot(334); ax(5) = subplot(335); ax(6) = subplot(336);
ax(7) = subplot(337); ax(8) = subplot(338); ax(9) = subplot(339);

for i=1:3
    set(ax(3*i),'tickdir','out')
    set(ax(3*i-1),'tickdir','out')
    set(ax(3*i-2),'tickdir','out')

    ylim(ax(3*i-2), [0 50])
    
    ylim(ax(3*i-1), [0 150])
    
    ylim(ax(3*i), [0 150])
end

slidingWindow(ax(1), cluster_bound_np, 1e1);
slidingWindow(ax(2), cluster_bound_tcr, 1e1);
slidingWindow(ax(3), cluster_phos_tcr, 1e1);
ylabel(ax(1), 'tau = 10')
title(ax(1), 'Bound NPs')
title(ax(2), 'Bound TCR')
title(ax(3), 'Phos TCR')

slidingWindow(ax(4), cluster_bound_np, 1e2);
slidingWindow(ax(5), cluster_bound_tcr, 1e2);
slidingWindow(ax(6), cluster_phos_tcr, 1e2);
ylabel(ax(4), 'tau = 100')

slidingWindow(ax(7), cluster_bound_np, 1e3);
slidingWindow(ax(8), cluster_bound_tcr, 1e3);
slidingWindow(ax(9), cluster_phos_tcr, 1e3);

ylabel(ax(7), 'tau = 1000')


%% Distributions of Adaptive Time Step

f4 = figure('Color','white');
ax(1) = subplot(231); ax(2) = subplot(232); ax(3) = subplot(233);
ax(4) = subplot(234); ax(5) = subplot(235); ax(6) = subplot(236);

for i=1:3 
    % Axis scale for Time Series
    set(ax(i),'tickdir','out')
    
    % Axis scales for Difference Series
    set(ax(i+3),'tickdir','out')
end

rho_vals = [1,100,31622];
for i = 1:3
    rho= rho_vals(i);
    
    load([file,'/Cluster_rho',num2str(rho),'.mat']);
    load([file,'/Uniform_rho',num2str(rho),'.mat']);
    
    TauDistribution(ax(i), cluster_bound_np, 'hist');
    xlabel(ax(i), 'Tau')
    title(ax(i), ['$\rho$ = ', num2str(rho)], 'Interpreter','Latex')

    TauDistribution(ax(i+3), homo_bound_np, 'hist');
    xlabel(ax(i+3), 'Tau')
    
end

%DifferenceSeries(ax(4), homo_bound_np);
%DifferenceSeries(ax(5), homo_bound_tcr);
%DifferenceSeries(ax(6), homo_phos_tcr);

%% Distributions of Adaptive Time Steps
function y = TauDistribution(ax, ts, plottype)

    for i = 1:size(ts,2)-4
        y = ts{i}(2,:);
        tau = y(2:end)-y(1:end-1);
        if strcmp(plottype, 'hist')
            histogram(ax, tau,'Normalization','pdf');
        elseif strcmp(plottype, 'series')
            plot(ax, [1:length(tau)], tau)
        else
            error(['Unrecognized plot type. Choose either `hist` or `series`.'])
        end
        hold(ax, 'on')
    end
end

%% Sliding Window Function

function f = slidingWindow(ax, ts, len_window)
    
    time_points = 1e5;
    end_time = 1e6;
    time_array = linspace(0,end_time,time_points);
    
    y = zeros(size(ts,2),time_points);
    f = zeros(1,time_points-len_window);
    
    for i=1:size(ts,2)
        y(i,:) = interp1(ts{i}(2,:), ts{i}(1,:), time_array);
    end
    
    len_final = time_points - len_window;
    
    j = randi([1 10],1,1);
    
    for i= 1:len_final
        window_start = i;
        window_end = i+len_window;
        f(i) = mean(y(:, window_start:window_end),'all');
    end
    hold(ax, 'on')
    plot(ax, time_array(1:end-len_window), f)
    xlabel(ax, 'Time')
end

%% Time Series Plotting Function
function y = TimeSeries(ax, ts)    
    for i=1:size(ts,2)
        hold(ax, 'on');
        plot(ax,ts{i}(2,:), ts{i}(1,:))
    end
    xlabel(ax, 'Time');
    ylim(ax, [0 300])
end

%% Difference Series Function
function y = DifferenceSeries(ax, ts)
    
    for i = 1:size(ts,2)
        y = ts{i}(1,:);
        t = ts{i}(2,2:end);
        Dy = y(2:end)-y(1:end-1);
        
        hold(ax,'on')
        plot(ax, t, Dy);
    end 
    
    xlabel(ax, 'Time');
    ylim(ax, [-100 100]);
end
%% Dickey-Fuller Plotting Function
function resultsDF = DFuller(ax, ts)
 
    h = zeros(size(ts,2), size(ts{1},2));
    
    for i=1:size(ts,2)
        y = ts{i}(1,:);
        t = ts{i}(2,:);
        total_t = t(end);
        
        final_time = linspace(0,total_t,1e3);
        
        Y = interp1(t, y, final_time, 'previous');
        
        for j = 111:length(Y)
            h(i,j) = adftest(Y(100:j));
        end
        
        %hold(ax,'on')
        %plot(ax,ts{i}(2,:), h(i,:));
    end
    
    imagesc(ax,h);
    ylabel(ax, 'Trial');
    xlabel(ax, 'Time');
    resultsDF = {'Bound NP', 'H', h};    
end

%%
%{

%% Autocorrelation Function

% Bound NPs
figure()
for i=1:size(cluster_bound_np,2)
    y = cluster_bound_np{i}(1,:);
    autocorr(y,NumLags=1e3-2);
    hold on
end

% Bound TCRs
figure()
for i=1:size(cluster_bound_tcr,2)
    y = cluster_bound_tcr{i}(1,:);
    autocorr(y, NumLags=1e3-2);
    hold on
end

% Phos TCRs
figure()
for i=1:size(cluster_phos_tcr,2)
    y = cluster_phos_tcr{i}(1,:);
    autocorr(y, NumLags=1e3-2);
    hold on
end

%% Time-Difference Series

% Bound NPs
figure()
for i=1:size(cluster_bound_np,2)
    y = cluster_bound_np{i}(1,:);
    yf = y(2:end) - y(1:end-1);
    plot([1:length(yf)], yf)
    hold on
end

% Bound TCRs
figure()
for i=1:size(cluster_bound_tcr,2)
    y = cluster_bound_tcr{i}(1,:);
    yf = y(2:end) - y(1:end-1);
    plot([1:length(yf)], yf)
    hold on
end

% Phos TCRs
figure()
for i=1:size(cluster_phos_tcr,2)
    y = cluster_phos_tcr{i}(1,:);
    yf = y(2:end) - y(1:end-1);
    plot([1:length(yf)], yf)
    hold on
end
%}