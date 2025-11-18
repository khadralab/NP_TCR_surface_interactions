function y = KPR_heatmap(ax1, valence, tpc, koff, kprtrue, negativefeedback)
    rho_vals = linspace(-4,3.5,16);
    x_vals = linspace(min(rho_vals), max(rho_vals), 1000);
    x_vals = 10.^x_vals;
    rho_vals = 10.^ rho_vals;
    koff_vals = [0.0001, 0.001, 0.01, 0.1];
    %koff_vals = [0.0001, 0.0003, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3];
    np_radius = 20;
    phos_rate = [0.0001, 0.001, 0.01, 0.1];
    phos_steps = [1, 2, 3, 4, 5];

    ec50_matrix = zeros(length(koff_vals), length(phos_rate));
    eMax_matrix = zeros(length(phos_rate), length(phos_steps));
    kpr_vals=zeros(length(koff_vals),length(x_vals));

    kp = [0.14, 0.07, 5]; ST=0.05; Cs=200;


    for i = 1:length(koff_vals)
        [mean_data, var_data, output_fit] = load_DR(koff_vals(i), [0.01, 0], rho_vals, tpc, 20, valence,1);
        fit_data = output_fit{1};
        y_vals = fit_data(x_vals);

        if ~kprtrue
            kpr_vals(i,:) = y_vals;
        else
            for j = 1:length(y_vals)
                kpr_vals(i,j) = negative_feedback_KPR(y_vals(j), koff_vals(i), kp, ST, Cs, negativefeedback);
                kpr_vals(i,j) = digital_activation(kpr_vals(i,j), 10);
            end
        end
    end
    %imagesc(ax1, x_vals, koff_vals, kpr_vals, [0 1]);
    [C,h] = contourf(ax1, x_vals, koff_vals, kpr_vals, [0 1 10 20 50 100],'ShowText','off','FaceAlpha',0.5,'LevelListMode','manual');
    %clabel(C,h,'LabelSpacing',72,'Margin',3,'Rotation',0);
    ax1.YDir = 'normal'; %ax2.Ydir = 'normal';
    %cb2 = colorbar(ax2, 'Ticks',[10^-4, 10^-3, 10^-2, 10^-1],'Location','eastoutside');
    %cb2.Label.String = 'EC50';
    set(ax1,'ColorScale','log');


end