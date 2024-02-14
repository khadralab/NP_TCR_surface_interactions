%% Dose-Response Curves
% Plots of results from DoseResponse.m for: (r,v, koff) = (20,5,0.05)
bound_type = {'phos', 'bound'};
file=['TauLeap'];

koff_vals = [0.01]; %, 0.02, 0.05, 0.1, 0.2, 0.5];
kd_vals = koff_vals ./ 0.1;
rho_vals = linspace(-2,5,15);           % Full range of NP-concentrations.
rho_vals = rho_vals(2:10);              % Select the x-axis concentrations to illustrate

f = figure('Color', 'white');
ax(1) = subplot(2,4,1:2); ax(3) = subplot(2,4,3); ax(5) = subplot(2,4,4);
ax(2) = subplot(2,4,5:6); ax(4) = subplot(2,4,7); ax(6) = subplot(2,4,8);

kd_lbls = arrayfun(@num2str, kd_vals, 'UniformOutput',0);

for i=1:2 
    % Axis scale for DR Curves
    set(ax(i),'xscale','log','tickdir','out')
    ylim(ax(i), [0 300])
    
    % Axis scales for Coop Plot
    set(ax(i+2),'xscale','log','yscale','log','tickdir','out')
    xlim(ax(i+2), [0.05 10]);
    ylim(ax(i+2), [0.1 10^4]);
    xticks(ax(i+2),[0.1, 1, 10]);
    
    set(ax(i+4),'xscale','log','yscale','log','tickdir','out')
    xlim(ax(i+4), [0.05 10]);
    ylim(ax(i+4), [10 10^3]);
    xticks(ax(i+4),[0.1, 1, 10]);
end

%Plotting

hols = autumn(length(koff_vals)+1);
hols = hols(1:length(koff_vals),:);
cols = winter(length(koff_vals));

scat1 = DoseResponsePlots(ax(1), kd_vals, rho_vals, 'clusters', bound_type{1}, hols, file);
scat2 = DoseResponsePlots(ax(2), kd_vals, rho_vals, 'uniform', bound_type{1}, cols, file);

CooperativityPlots(ax(3), kd_vals, rho_vals, 'clusters', bound_type{1}, hols, file)
CooperativityPlots(ax(4), kd_vals, rho_vals, 'uniform', bound_type{1}, cols, file)

EMaxPlots(ax(5), kd_vals, rho_vals, 'clusters', bound_type{1}, hols, file)
EMaxPlots(ax(6), kd_vals, rho_vals, 'uniform', bound_type{1}, cols, file)

% Legend
l1 = legend(ax(1), scat1, kd_lbls,'location','northwest','NumColumns',2);
l2 = legend(ax(2), scat2, kd_lbls,'location','northwest','NumColumns',2);
title(l1, 'Clusters $K_D$');
title(l2, 'Uniform $K_D$');

set(findall(f,'-property','FontSize'),'FontSize',16)
set(findall(f,'-property','Interpreter'),'Interpreter','Latex')

%% Hill Coefficient Figure

f2 = figure('Color','white');
ax(1) = subplot(121); ax(2) = subplot(122);

for i=1:2
    set(ax(i),'xscale','log','tickdir','out')
    xlim(ax(i), [0.05 10]);
    ylim(ax(i), [0 2]);
    xticks(ax(i),[0.1, 1, 10]);
end

%Plotting

hols = autumn(length(koff_vals)+1);
hols = hols(1:length(koff_vals),:);
cols = winter(length(koff_vals));

HillCoeffPlots(ax(1), kd_vals, rho_vals, 'clusters', bound_type{1}, hols, file)
HillCoeffPlots(ax(2), kd_vals, rho_vals, 'uniform', bound_type{1}, cols, file)

set(findall(f2,'-property','FontSize'),'FontSize',16)
set(findall(f2,'-property','Interpreter'),'Interpreter','Latex')

%% Local Functions
% Cooperativity
function y = CooperativityPlots(ax, kd_vals, rho_vals, surface_type, bound_type, colors, file)
    x_vals = linspace(min(rho_vals), max(rho_vals), 1000);
    rho_vals = 10.^rho_vals;
    x_vals = 10.^x_vals;
    
    bound50 = []; phos50 = [];
    
    i=1;
    
    for koff = kd_vals ./ 10
        [mean_bound, mean_phos, var_bound, var_phos, fit_bound, fit_phos] = load_DR(koff, rho_vals, surface_type, file);
        bound50 = [bound50,fit_bound.ec50];
        phos50 = [phos50, fit_phos.ec50];
    end
    
    if strcmp(bound_type,'phos')
        ec50 = phos50;
    else
        ec50 = bound50;
    end
    
    x = log10(kd_vals); y = log10(ec50)';
    X = [ones(length(x),1),x'];

    b = X \ y;
    
    hold(ax,'on');
    scatter(ax, kd_vals, ec50, [], colors, 'filled');
    plot(ax, kd_vals, 10^b(1)*kd_vals.^b(2));
    title(ax, ['Intercept = ',num2str(b(1)),', Slope = ',num2str(b(2))]);
    ylabel(ax,'EC50');
    xlabel(ax,'$K_D$');
    grid(ax,'on');
    
end

%% EMAX Plots

function y = EMaxPlots(ax, kd_vals, rho_vals, surface_type, bound_type, colors, file)
    x_vals = linspace(min(rho_vals), max(rho_vals), 1000);
    rho_vals = 10.^rho_vals;
    x_vals = 10.^x_vals;
    
    boundMax = []; phosMax = [];
    
    i=1;
    
    for koff = kd_vals ./ 10
        [mean_bound, mean_phos, var_bound, var_phos, fit_bound, fit_phos] = load_DR(koff, rho_vals, surface_type, file);
        boundMax = [boundMax,fit_bound.eMax];
        phosMax = [phosMax, fit_phos.eMax];
    end
    
    if strcmp(bound_type,'phos')
        eMax = phosMax;
    else
        eMax = boundMax;
    end
    
    x = log10(kd_vals); y = log10(eMax)';
    X = [ones(length(x),1),x'];

    b = X \ y;
    
    hold(ax,'on');
    scatter(ax, kd_vals, eMax, [], colors, 'filled');
    plot(ax, kd_vals, 10^b(1)*kd_vals.^b(2));
    title(ax, ['Intercept = ',num2str(b(1)),', Slope = ',num2str(b(2))]);
    ylabel(ax,'EMax');
    xlabel(ax,'$K_D$');
    grid(ax,'on');
    
end

%% Hill Coefficient Plots

function y = HillCoeffPlots(ax, kd_vals, rho_vals, surface_type, bound_type, colors, file)
    x_vals = linspace(min(rho_vals), max(rho_vals), 1000);
    rho_vals = 10.^rho_vals;
    x_vals = 10.^x_vals;
    
    boundHill = []; phosHill = [];
    
    i=1;
    
    for koff = kd_vals ./ 10
        [mean_bound, mean_phos, var_bound, var_phos, fit_bound, fit_phos] = load_DR(koff, rho_vals, surface_type, file);
        boundHill = [boundHill,fit_bound.n];
        phosHill = [phosHill, fit_phos.n];
    end
    
    if strcmp(bound_type,'phos')
        hill_cof = phosHill;
    else
        hill_cof = boundHill;
    end
    
    x = log10(kd_vals); y = hill_cof';
    X = [ones(length(x),1),x'];

    b = X \ y;
    
    hold(ax,'on');
    scatter(ax, kd_vals, hill_cof, [], colors, 'filled');
    plot(ax, kd_vals, log10(10.^b(1)*kd_vals.^b(2)));
    title(ax, ['Intercept = ',num2str(b(1)),', Slope = ',num2str(b(2))]);
    ylabel(ax,'Hill Coefficient');
    xlabel(ax,'$K_D$');
    grid(ax,'on');
end

%% Dose Response Plots

function [scat] = DoseResponsePlots(ax, kd_vals, rho_vals, surface_type, bound_type, colors, file)
    x_vals = linspace(min(rho_vals), max(rho_vals), 1000);
    rho_vals = 10.^rho_vals;
    x_vals = 10.^x_vals;
    
    i=1;

    for koff = kd_vals ./ 10
        [mean_bound, mean_phos, var_bound, var_phos, fit_bound, fit_phos] = load_DR(koff, rho_vals, surface_type, file);
        std_bound = sqrt(var_bound); std_phos = sqrt(var_phos);
       
        % Plot mean Phosphorylated dose-response, Hill fit and shade standard deviation
        if strcmp(bound_type,'phos')
            hold(ax, 'on');
            scat(i) = scatter(ax,rho_vals, mean_phos,[], colors(i,:),'filled', 'DisplayName',['Cluster ',num2str(koff)]);
            plot(ax,x_vals, fit_phos(x_vals), 'Color', colors(i,:));
            shade(ax, rho_vals, mean_phos+std_phos, 'w', rho_vals, max(0,mean_phos-std_phos),'w', 'FillType',[1 2;2 1],'FillColor', colors(i,:),'FillAlpha',0.1, 'HandleVisibility', 'off');
            grid(ax, 'on');
            ylabel(ax,'Phos. TCRs');
            xlabel(ax,'NP Concentration');
        else
            hold(ax, 'on');
            scat(i) = scatter(ax,rho_vals, mean_bound,[], colors(i,:),'filled', 'DisplayName',['Cluster ',num2str(koff)]);
            plot(ax,x_vals, fit_bound(x_vals), 'Color', colors(i,:));
            shade(ax, rho_vals, mean_bound+std_bound, 'w', rho_vals, max(0,mean_bound-std_bound),'w', 'FillType',[1 2;2 1],'FillColor', colors(i,:),'FillAlpha',0.1, 'HandleVisibility', 'off');
            grid(ax, 'on');
            ylabel(ax,'Bound TCRs');
            xlabel(ax,'NP Concentration');
        end
        
        i=i+1;
    end
end

%% Load Data Function

function [mean_bound, mean_phos, var_bound, var_phos, fit_bound, fit_phos] = load_DR(koff, rho_vals, surface_type, file)
    num_sims = 10;
    mean_bound = []; var_bound = [];
    mean_phos = []; var_phos = [];
    
    f = ['LongSims/',file,'/'];
    
    mean_bt = zeros(1,num_sims); mean_pt = zeros(1,num_sims);
        
    if strcmp(surface_type,'clusters')
        for rho = rho_vals
            load([f,'koff',num2str(koff*100),'/Cluster_rho',num2str(floor(rho*100)),'.mat'])
            
            for i=1:num_sims
                bt = cluster_bound_tcr{i}(1,end-2000:end); pt = cluster_phos_tcr{i}(1,end-2000:end);
                bt = bt(~isnan(bt)); pt = pt(~isnan(pt));
                mean_bt(i) = mean(bt); mean_pt(i) = mean(pt);
            end

            mean_bound = [mean_bound, mean(mean_bt)];
            var_bound = [var_bound, var(mean_bt)];
            mean_phos = [mean_phos, mean(mean_pt)];
            var_phos = [var_phos, var(mean_pt)];
        end
        
    else
        for rho = rho_vals
            load([f,'koff',num2str(koff*100),'/Uniform_rho',num2str(floor(rho*100)),'.mat'])

            for i=1:num_sims
                bt = homo_bound_tcr{i}(1,end-2000:end); pt = homo_phos_tcr{i}(1,end-2000:end);
                bt = bt(~isnan(bt)); pt = pt(~isnan(pt));
                mean_bt(i) = mean(bt); mean_pt(i) = mean(pt);
            end

            mean_bound = [mean_bound, mean(mean_bt)];
            var_bound = [var_bound, var(mean_bt)];
            mean_phos = [mean_phos, mean(mean_pt)];
            var_phos = [var_phos, var(mean_pt)];

        end
    end
    
    % Fit Hill Function to dose response

    fo = fitoptions('Method','NonlinearLeastSquares',...
                   'Lower',[0,0],...
                   'Upper',[Inf,max(rho_vals)],...
                   'StartPoint',[1 1]);
    ft = fittype('eMax*x^n / (x^n + ec50^n)','problem','n','options',fo);

    [fit_bound,gof] = fit(rho_vals',mean_bound',ft,'problem',1);
    [fit_phos,gof] = fit(rho_vals',mean_phos',ft,'problem',1);
    
end

%% Fitting Function

function [meanC, meanU, varC, varU, c_ft, u_ft] = Hill_fitting(koff, rho_vals, surface_type)

    [meanC, meanU, varC, varU] = DR_fitting(koff, rho_vals, surface_type);
    

    % Fit Hill Function to dose response

    fo = fitoptions('Method','NonlinearLeastSquares',...
                   'Lower',[0,0],...
                   'Upper',[Inf,max(rho_vals)],...
                   'StartPoint',[1 1]);
    ft = fittype('eMax*x^n / (x^n + ec50^n)','problem','n','options',fo);

    [c_ft,c_gof] = fit(rho_vals',meanC(2,:)',ft,'problem',0.5);

    [u_ft,u_gof] = fit(rho_vals',meanU(2,:)',ft,'problem',1);

end