%% Distributions of Bound TCRs
clear all

file = ['LongSims/TauLeap/v2/koff1e3/'];
surface = ['uniform'];

rho_vals = linspace(-2,5,15);
rho_vals = 10.^rho_vals;

rho = rho_vals(10);
total_t = 1e6;
wind = 1e3;

M = Distribution_Movie(rho, total_t, wind, file, surface);

movie(gcf, M, 1, 10);


%%

function [M] = Distribution_Movie(rho, total_t, wind, file, surface)
    n_frames = total_t / wind;
    final_time = linspace(0,total_t,n_frames);
    bound_tcr = load_data(file, rho, n_frames, wind, surface);
    
    f = figure('Color', 'white');
    f.Visible = 'off';
    ax1 = subplot(2,2,1:2);
    ax2 = subplot(223);
    ax3 = subplot(224);

    xlim(ax1, [0 300])
    ylim(ax1, [0 0.3])
    
    edges = linspace(0,300,101);
    
    ax1.NextPlot = 'replaceChildren';
    M(n_frames) = struct('cdata',[],'colormap',[]);
    
    
    for i = 1:n_frames
        [alpha(i), beta(i)] = Distributions(ax1, bound_tcr{i}, edges);
        plot(ax2, [1:i], alpha.*beta)
        xlabel(ax2, 'Time')
        ylabel(ax2, 'Mean')
        
        plot(ax3, [1:i], alpha.*beta.^2)
        xlabel(ax3, 'Time')
        ylabel(ax3, 'Variance')
        drawnow;
        M(i) = getframe(gcf);
    end
    f.Visible = 'on';
end

%%
function [alpha, beta] = Distributions(ax, bound_tcr, edges)
    histogram(ax, bound_tcr, 'BinEdges', edges, 'Normalization', 'pdf');
    hold(ax, 'on')
    x_vals = linspace(0,300,301);
    pd = fitdist(bound_tcr','gamma');
    alpha = pd.a;
    beta = pd.b;
    y = pdf(pd, x_vals);
    plot(ax,x_vals, y)
    hold(ax,'off')
    xlabel(ax,'Bound TCRs')
    ylabel(ax,'Probability Density')
    drawnow;
end
%%
function [bound_tcr] = load_data(file, rho, n_frames, wind, surface)
    bound_tcr{1, n_frames} = {};
    
    if strcmp(surface,'cluster')
        load([file, 'Cluster_rho',num2str(floor(rho*100)),'.mat']);

        for i =1:length(cluster_bound_tcr)
            bt = cluster_bound_tcr{i}(1,:);
            t = cluster_bound_tcr{i}(2,:);
            
            for j=1:n_frames
                c = (t>= wind*(j-1)) & (t < wind*j);
                
                bound_tcr{j} = [bound_tcr{j}, bt(c)];
            end
        end     
            
    else
        load([file, 'Uniform_rho',num2str(floor(rho*100)),'.mat']);
        for i =1:length(homo_bound_tcr)
            bt = homo_bound_tcr{i}(1,:);
            t = homo_bound_tcr{i}(2,:);
            
            for j=1:n_frames
                c = (t>= wind*(j-1)) & (t < wind*j);
                
                bound_tcr{j} = [bound_tcr{j}, bt(c)];
            end
        end  
    end
    
    bound_tcr{n_frames} = cell2mat(bound_tcr{n_frames});
end