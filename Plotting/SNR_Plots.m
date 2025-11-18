%% SNR Plots
clear all

bound_type = {'bound', 'phos', 'np'};
btype = 1;
tpc = 0;

rho_vals = linspace(-4,3.5,16);
rho_vals = 10.^rho_vals;

rho = [rho_vals(11), rho_vals(11), rho_vals(11), rho_vals(11)];
koff_vals = [0.0001, 0.001, 0.01, 0.1];
valence = [5,5,5,5];

kd_lbls = arrayfun(@num2str, koff_vals ./ 0.1, 'UniformOutput',0);

cols = cool(length(koff_vals));
hols = hot(length(koff_vals)+4);
gols = gray(length(koff_vals)+1);
aols = abyss(length(koff_vals));

f1 = figure('Color', 'white','Units','centimeters','Position',[10 10 37 12]);
f2 = figure('Color', 'white','Units','centimeters','Position',[10 10 37 12]);
f3 = figure('Color', 'white','Position',[10 100 1500 375]);

left = [0.1, 0.37, 0.07, 0.55, 0.07, 0.55, 0.07, 0.55];
bottom = 0.4;
width = 0.21;
height =0.5;

for i=1:2
    ax(i) = axes(f1, 'Units','normalized','Position',[left(i) bottom width height]);
end
for i=3:4
    ax(i) = axes(f2, 'Units','normalized','Position',[left(i) 0.3 0.25 height]);
end
for i=5:6
    ax(i) = axes(f3, 'Units','normalized','Position',[left(i) bottom width height]);
end

ax(7) = axes(f1, 'Units','normalized','Position', [0.64 bottom width height]);
hold(ax(:),'on')

% Specificity (Kd)
[h1,a1,v1] = TCR_Distributions(ax(2), koff_vals, rho, valence, 3, bound_type{1},aols);
[h2,a2,v2] = TCR_Distributions(ax(1), koff_vals, rho, valence, 0, bound_type{1},aols);
[h7,a7,v7] = TCR_Distributions(ax(7), koff_vals, rho, valence, 20, bound_type{1},aols);

%text(ax(2), a1+[5,-22,-22,5], [0.06 0.035 0.03 0.04], arrayfun(@num2str, floor(a1),'UniformOutput',0))
%text(ax(1), a2+[3,7,7,3], [0.04 0.03 0.06 0.12], arrayfun(@num2str, floor(a2),'UniformOutput',0))
%text(ax(7), a7+[7,-22,10,3], [0.05 0.03 0.02 0.06], arrayfun(@num2str, floor(a7),'UniformOutput',0))

% Axis labels and format
set([ax(1:4) ax(7)],'GridLineWidth',0.5,'Linewidth',2)
xlim(ax(1:2), [0, 300]);
xlim(ax(7), [0, 300]);
yticks(ax(2), [0, 0.02, 0.04, 0.06])
yticks(ax(7), [0, 0.02, 0.04, 0.06])
xlabel([ax(1) ax(7)], '')

% Valence sensitivity
rho = [rho_vals(15), rho_vals(15), rho_vals(15), rho_vals(15)];
koff_vals = [0.001, 0.001, 0.001, 0.001];
valence = [1,2,5,10];
val_lbl = arrayfun(@num2str,valence,'UniformOutput',0);

[h3,a3,v3] = TCR_Distributions(ax(3), koff_vals, rho, valence, 20, bound_type{1},hols);
[h4,a4,v4] = TCR_Distributions(ax(4), koff_vals, rho, valence, 0, bound_type{1},gols);

%text(ax(3), a3+[5,5,-20,5], [0.12 0.04 0.04 0.04], arrayfun(@num2str, floor(a3),'UniformOutput',0))
%text(ax(4), a4+[3,3,3,3], [0.12 0.05 0.05 0.05], arrayfun(@num2str, floor(a4),'UniformOutput',0))

ylabel(ax(4), 'Probability density');

% Dose sensitivity
rho = [rho_vals(2), rho_vals(6), rho_vals(10), rho_vals(14)];
koff_vals = [0.0001, 0.0001, 0.0001, 0.0001];
valence = [5,5,5,5];

h5 = TCR_Distributions(ax(5), koff_vals, rho, valence, 20, bound_type{1},hols);
h6 = TCR_Distributions(ax(6), koff_vals, rho, valence, 0, bound_type{1},gols);

for i=1:3
    ylabel(ax(2*i-1),'Probability density');
    %xlim(ax(i*2-1), [0 300]);
    %xlim(ax(i*2), [0 300]);
end

set(ax(:),'tickdir','out','box','off');
grid(ax(:),'on');

% Legend
%pos1 = get(ax(5),'Position');
%pos2 = get(ax(6),'Position');

l1 = legend(ax(1), h1, kd_lbls,'Position',[0.88, 0.6, 0.1, 0.1],'NumColumns',1);
l2 = legend(ax(3), h3, val_lbl,'Position',[0.34, 0.47, 0.1, 0.1],'NumColumns',1);
l3 = legend(ax(4), h4, val_lbl,'Position',[0.82, 0.47, 0.1, 0.1],'NumColumns',1);
l4 = legend(ax(5), h5, arrayfun(@num2str,log10(rho),'UniformOutput',0),'Position',[0.43, 0.4, 0.1, 0.1],'NumColumns',2);
title(l1, '$K_D$');
title(l2, 'NP Valence');
title(l3, 'NP Valence');
title(l4, 'log10(Concentration)');


% Set fontsize
set(findall(f1,'-property','FontSize'),'FontSize',20)
set(findall(f1,'-property','Interpreter'),'Interpreter','Latex')
set(findall(f2,'-property','FontSize'),'FontSize',20)
set(findall(f2,'-property','Interpreter'),'Interpreter','Latex')
set(findall(f3,'-property','FontSize'),'FontSize',20)
set(findall(f3,'-property','Interpreter'),'Interpreter','Latex')

title(ax(7),' 20 TPC','Color','r','FontSize',20);
title(ax(1), ' 1 TPC (Uniform)','Color','k','FontSize',20);
title(ax(2), '3 TPC','Color','b','FontSize',20);

% Save
set(f2, 'Units', 'inches')
pos = get(f2, 'Position');
set(f2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(f2, '../Figures/Manuscript/Final-Figures/Valence-Distributions','-dpdf','-r0','-fillpage');

%%

f = figure('Color', 'white');

ax(1) = subplot(2,2,1:2); ax(2) = subplot(2,2,3:4);

valence = 5;
tpc = [3];
cols = winter(length(tpc));

rho_vals = linspace(-4,3.5,16);
rho_vals = 10.^rho_vals;
koff_vals = [0.3, 0.1, 0.03, 0.01, 0.003, 0.001, 0.0003, 0.0001];
%koff_vals = [0.1, 0.01, 0.001, 0.0001];

for i=1:length(tpc)
    Plot_dprime(ax(1), [koff_vals(1), koff_vals(1)], rho_vals, [valence, valence], tpc(i), bound_type, cols, true)
    hold(ax(1),'on');

    Plot_dprime(ax(2), koff_vals, [rho_vals(8), rho_vals(8)], [valence, valence], tpc(i), bound_type, cols, false)
    hold(ax(2),'on')
end

l1 = legend(ax(1),{'Uniform', '3', '5', '10', '20'});
title(l1,'TCRs per Cluster')

%% Local Functions
%% D-Prime

function [F] = Plot_dprime(ax, koff_vals, rho_vals, valence, tpc, bound_type, cols, arg)
    if arg
        x(:,1) = rho_vals(1:end-1)';
        x(:,2) = rho_vals(2:end)';
        
        xlbl = 'Rho';
        
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

    semilogx(ax,x(:,1),dp)
    ylabel(ax,' D-Prime')
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

%% Receiver Operating Characteristic

function [F] = ROC(ax, koff_vals, rho, valence, tpc, bound_type, cols)

    data_train = []; data_test = [];
    label_train = []; label_test = [];
    
    for i = 1:2
        
        koff= koff_vals(i);
        
        [tcr_dist] = load_DR(koff, rho(i), valence(i), tpc, bound_type);
        
        start = floor(0.99*length(tcr_dist));
        start2 = floor(0.98*length(tcr_dist));
        
        x = reshape(tcr_dist(:,start:end),[],1);
        x2 = reshape(tcr_dist(:,start2:start),[],1);
        
        y = i * ones(size(x));
        y2 = i * ones(size(x2));
        
        data_train = [data_train; x];
        label_train = [label_train; y];
        
        data_test = [data_test; x2];
        label_test = [label_test; y2];
    end
    
    idx = randperm(length(data_train));
    data_train = data_train(idx); label_train = label_train(idx);
    
    idx = randperm(length(data_test));
    data_test = data_test(idx); label_test = label_test(idx);
    
    mdl = fitglm(data_train,label_train,'Distribution','poisson','Link','log');
    
    predictions = predict(mdl, data_test);
    [X,Y,T,AUC] = perfcurve(label_test,predictions, 2);
    
    plot(ax,X,Y)
    xlabel(ax,'False positive rate') 
    ylabel(ax,'True positive rate')
    title(ax, AUC);
end

%% Plotting Histograms

function [H, average, variance] = TCR_Distributions(ax, koff_vals, rho_vals, valence, tpc, bound_type, cols)
    hold(ax,'on');
    
    for i=1:length(koff_vals)
        koff = koff_vals(i); rho = rho_vals(i); v = valence(i);
        
        [tcr_dist] = load_DR(koff, rho, v, tpc, bound_type);
        
        start = floor(0.75*length(tcr_dist));
        
        average(i) = mean(tcr_dist(:,start:end),'all');
        variance(i) = var(tcr_dist(:,start:end),0,'all');
    
        H(i) = histogram(ax, tcr_dist(:,start:end),'FaceAlpha',0.8,'EdgeAlpha',0.2,'FaceColor',cols(i,:),'Normalization','probability','BinMethod','integers');
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

function [tcr_dist] = load_DR(koff, rho, v, tpc, bound_type)

    file = ['LongSims/r20/',num2str(tpc),'TPC/v',num2str(v),'/'];
    
    if tpc == 0 
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

%% Bootstrap Method

function [snr, mean_snr, var_snr] = Bootstrap_SNR(koff, rho, valence, tpc, bound_type)
    tcr_dist = load_DR(koff, rho, valence, tpc, bound_type);
    
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