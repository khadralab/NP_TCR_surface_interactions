load(['MixtureResults/koff1-5/Cluster_rho1-50.mat'])

figure('Color', 'white')
ax(1) = subplot(221);
ax(2) = subplot(222);
ax(3) = subplot(223);
ax(4) = subplot(224);

for i=1:10
    plot(ax(1), cluster_bound_np1{i}(2,:), cluster_bound_np1{i}(1,:))
    hold(ax(1),'on')
    plot(ax(2), cluster_bound_np2{i}(2,:), cluster_bound_np2{i}(1,:))
    hold(ax(2),'on')
    plot(ax(3), cluster_bound_tcr{i}(2,:), cluster_bound_tcr{i}(1,:))
    hold(ax(3),'on')
    plot(ax(4), cluster_phos_tcr{i}(2,:), cluster_phos_tcr{i}(1,:))
    hold(ax(4),'on')
end

title(ax(1), 'NP Type 1')
title(ax(2), 'NP Type 2')
title(ax(3), 'Bound TCRs')
title(ax(4), 'Phos TCRs')



%%

f = figure('Color', 'white');
ax(1) = subplot(221);
ax(2) = subplot(222);
ax(3) = subplot(223);
ax(4) = subplot(224);

for i=1:10
    plot(ax(1), homo_bound_np1{i}(2,:), homo_bound_np1{i}(1,:))
    hold(ax(1),'on')
    plot(ax(2), homo_bound_np2{i}(2,:), homo_bound_np2{i}(1,:))
    hold(ax(2),'on')
    plot(ax(3), homo_bound_tcr{i}(2,:), homo_bound_tcr{i}(1,:))
    hold(ax(3),'on')
    plot(ax(4), homo_phos_tcr{i}(2,:), homo_phos_tcr{i}(1,:))
    hold(ax(4),'on')
end

title(ax(1), 'NP Type 1')
title(ax(2), 'NP Type 2')
title(ax(3), 'Bound TCRs')
title(ax(4), 'Phos TCRs')
%%
clear all
file = ['MixtureResults/koff1-5/'];

f = figure('Color', 'white');
ax1 = subplot(121);
ax2 = subplot(122);

plotMix(ax1, 1, linspace(0,1,11), file, 'cluster');
title(ax1, ['Cluster Surface'])

plotMix(ax2, 1, linspace(0,1,11), file, 'uniform');
title(ax2, ['Uniform Surface'])

set(findall(f,'-property','FontSize'),'FontSize',24)
set(findall(f,'-property','Interpreter'),'Interpreter','Latex')
%%

function y = plotMix(ax, rho, mix_vals, file, surface)
    [tcr, np1, np2] = load_data(rho, mix_vals, file, surface);
    
    % Right axis
    yyaxis(ax,'right')
    plot(ax, mix_vals, tcr,'r-', 'Linewidth',2)
    ylabel(ax, 'Bound TCR')
    hold(ax, 'on')
    ylim(ax, [0 5])
    
    % Left Axis
    yyaxis(ax,'left')
    plot(ax, mix_vals, np1,'b-','Linewidth',2)
    plot(ax, mix_vals, np2,'g-','Linewidth',2)
    ylabel(ax, 'Bound NPs')
    xlabel(ax, 'NP Ratio')
    ylim(ax, [0 5])
    
    %xline(ax,0.3, '--k', 'Linewidth',2)
    
    set(ax,'XMinorTick','on','tickdir','out')
    grid(ax,'on')
    grid(ax, 'minor')
    l = legend(ax, {'$k_{off}=0.01$', '$k_{off}=0.05$','','TCR'},'location','northwest');
    
end

function [bound_tcr, bound_np1, bound_np2] = load_data(rho, mix_vals, file, surface)
    bound_np1 = [];
    bound_np2 = [];
    bound_tcr = [];
    
    if strcmp(surface,'cluster')
        fname = [file,'Cluster_rho'];
        
        for i = 1:length(mix_vals)
            mixRatio = mix_vals(i);
            load([fname,num2str(rho),'-',num2str(mixRatio*100),'.mat']);
            
            for j=1:length(cluster_bound_np1)
                bt = cluster_bound_tcr{j}(1,end-2000:end);
                np1 = cluster_bound_np1{j}(1,end-2000:end);
                np2 = cluster_bound_np2{j}(1,end-2000:end);
                
                bt = bt(~isnan(bt)); np1 = np1(~isnan(np1)); np2 = np2(~isnan(np2));
                mean_bt(j) = mean(bt); mean_np1(j) = mean(np1); mean_np2(j) = mean(np2);
            end
            
            bound_tcr = [bound_tcr, mean(mean_bt)];
            bound_np1 = [bound_np1, mean(mean_np1)];
            bound_np2 = [bound_np2, mean(mean_np2)];

        end
    else
        fname = [file,'Uniform_rho'];
        
        for i = 1:length(mix_vals)
            mixRatio = mix_vals(i);
            load([fname,num2str(rho),'-',num2str(mixRatio*100),'.mat']);
            
            for j=1:length(homo_bound_np1)
                bt = homo_bound_tcr{j}(1,end-2000:end);
                np1 = homo_bound_np1{j}(1,end-2000:end);
                np2 = homo_bound_np2{j}(1,end-2000:end);
                
                bt = bt(~isnan(bt)); np1 = np1(~isnan(np1)); np2 = np2(~isnan(np2));
                mean_bt(j) = mean(bt); mean_np1(j) = mean(np1); mean_np2(j) = mean(np2);
            end
            
            bound_tcr = [bound_tcr, mean(mean_bt)];
            bound_np1 = [bound_np1, mean(mean_np1)];
            bound_np2 = [bound_np2, mean(mean_np2)];

        end
    end

    
    
end



