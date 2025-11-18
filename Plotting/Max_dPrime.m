% Max d-prime Plots

clear all
addpath('Plotting/Functions')

bound_type = {'phos', 'bound', 'np', 'covered','fracphos', 'xPMHC', 'rescaled'};
btype = 1;                          % 1 for phos; 2 for bound; 3 for NPs
tcrs_per_cluster = [0, 3, 5, 10];
valence = 5;

% Kd and NP rho
koff_vals = [0.0001, 0.0003, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3];
kd_vals = koff_vals ./ 0.1;
rho_vals = linspace(-4,3.5,16);

f = figure('Color','white','Position',[10 100 1500 800]);
ax(1) = subplot(321); ax(2) = subplot(322); ax(3) = subplot(323); ax(4) = subplot(324); ax(5) = subplot(325); ax(6) = subplot(326);
a1_pos = get(ax(1), 'Position');
a2_pos = get(ax(2), 'Position');

hols = autumn(length(koff_vals)+1);
hols = hols(1:length(koff_vals),:);
cols = winter(length(rho_vals));

kd_lbls = arrayfun(@num2str, kd_vals, 'UniformOutput',0);
xlbl = 10.^(-3:3);
xlabels = arrayfun(@num2str, xlbl, 'UniformOutput',0);

% Axis for D-Prime Plot
set(ax(1:6),'xscale','log','tickdir','out')
set(ax([4,6]), 'yscale','log')

rho_vals = 10.^rho_vals;

for TPC = tcrs_per_cluster
    clear max_dp
    for i=1:length(koff_vals)
        hold(ax(1),'on')
        [dp, x] = Plot_dprime(ax(1), [koff_vals(i), koff_vals(i)], rho_vals, [valence, valence], TPC, bound_type{btype}, hols(i,:), true);
        
        max_dp(i) = max(dp);
        conc_at_max(i) = x(dp == max_dp(i));
    end
    
    hold(ax(3:6), 'on')
    plot(ax(3), kd_vals, max_dp, ':x', 'Linewidth', 1.5, 'DisplayName', num2str(TPC))
    plot(ax(4), kd_vals, conc_at_max, ':x', 'Linewidth', 1.5, 'DisplayName', num2str(TPC))

    clear max_dp
    for i=1:length(rho_vals)
        hold(ax(2),'on')
        [dp, x] = Plot_dprime(ax(2), koff_vals, [rho_vals(i), rho_vals(i)], [valence, valence], TPC, bound_type{btype}, cols(i,:), false);
        max_dp(i) = max(dp);
        kd_at_max(i) = x(dp == max_dp(i));
    end

    plot(ax(5), rho_vals, max_dp, ':x', 'Linewidth', 1.5, 'DisplayName', num2str(TPC))
    plot(ax(6), rho_vals, kd_at_max, ':x', 'Linewidth', 1.5, 'DisplayName', num2str(TPC))
    
end

xlabel(ax(3:4), '$K_D$')
ylabel(ax(3), 'Max d-prime')
ylabel(ax(4), 'NP concentration at maximum d-prime')

xlabel(ax(5:6), 'NP concentration')
ylabel(ax(5), 'Max d-prime')
ylabel(ax(6), '$K_D$ at maximum d-prime')

l1 = legend(ax(3)); l2 = legend(ax(4));

set(findall(f,'-property','FontSize'),'FontSize',14)
set(findall(f,'-property','Interpreter'),'Interpreter','Latex')