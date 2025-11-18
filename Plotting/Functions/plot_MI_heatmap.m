function y = plot_MI_heatmap(ax,TPC,mi_variable,savetrue,type)

    %fpath = ["/NP-Revival/SurfaceMovie/"];
    %fpath = ["/media/louis/Untitled/NP-Revival/SurfaceMovie/"];
    fpath = pwd+"/Mutual Information/";
    if strcmp(type, "norm")
        fend = ["_nmi.mat"];
    elseif strcmp(type, "entropy")
        fend = ["_entropy.mat"];
    else
        fend = ["_mi.mat"];
    end

    valence_vals = [1,2,3,4,5];
    rho_vals = linspace(-4,3.5,16);
    rho_vals = 10.^ rho_vals;
    koff_vals = [0.0001, 0.001, 0.01, 0.1];
    surfaces = [1,3,5,10,20];

    t = find(surfaces == TPC);

    if strcmp(mi_variable,"kd")
        x_array = rho_vals;
        xlbl = 'NP Concentration';
        xtf = '%1.f';
        xt = log10([10^-3, 10^-1, 10^1, 10^3]);%log10(rho_vals(1:end));
        y_array = 10.^valence_vals;
        ylbl = 'Valence';
        yt = valence_vals;
        ytf='%1.f';
        load(fpath+"kd"+fend);
        data = mi_cell{1,t};
        data = data(1:length(y_array),1:length(x_array));

    elseif strcmp(mi_variable, "valence")
        x_array = rho_vals;
        xlbl = 'NP Concentration';
        xt = log10([10^-3, 10^-1, 10^1, 10^3]); %log10(rho_vals(1:end));
        y_array = koff_vals * 10;
        ylbl = '$K_D$';
        yt = log10(koff_vals * 10);
        load(fpath+"valence"+fend);
        ytf='10^{%1.f}';
        xtf='10^{%1.f}';
        data = mi_cell{1,t}';
        data = data(1:length(y_array),1:length(x_array));

    elseif strcmp(mi_variable, "concentration")
        x_array = koff_vals*10;
        xlbl = '$K_D$';
        xt = log10(koff_vals*10);
        xtf = '10^{%1.f}';
        y_array = 10.^valence_vals;
        yt = valence_vals;
        ytf = '%1.f';
        ylbl = 'Valence';
        load(fpath+"rho"+fend);
        data = mi_cell{1,t};
        data = data(1:length(y_array),1:length(x_array));

    else
        error("Variable must be either 'kd', 'valence' or 'concentration'. ")
    end

    channel_cap = max(data,[],'all');
    channel_sum = mean(data, "all");
    
    colormap(jet);
    
    %levels = [0.,0.25,0.5,0.75,1.,1.25,1.5,1.75,2.,2.25];
    %contourf(ax,log10(x_array), log10(y_array), data,'ShowText','on','FaceAlpha',0.55,'LevelList', levels)
    
    if strcmp(type, "norm")
        imagesc(ax, log10(x_array), log10(y_array), data, [0 0.27])
        ax.YDir = 'normal';
        if TPC == 20
            cb = colorbar(ax, 'Ticks',[0,0.05,0.1,0.15,0.2,0.25]);
            cb.Label.String = 'Mutual Information (bits)';
        end
    elseif strcmp(type, "entropy")
        imagesc(ax, log10(x_array), log10(y_array), data, [0 8])
        ax.YDir = 'normal';
        if TPC == 20
            cb = colorbar(ax, 'Ticks',[0 2 4 6 8]);
            cb.Label.String = 'Entropy';
        end
    else
        imagesc(ax, log10(x_array), log10(y_array), data, [0 2.5])
        ax.YDir = 'normal';
        if TPC == 20
            cb = colorbar(ax, 'Ticks',[0,0.5,1,1.5,2,2.5],'Location','eastoutside');
            cb.Label.String = 'Mutual Information (bits)';
        end
    end

    grid(ax,'off');
    %ax.YMinorGrid='on';
    %grid(ax,'minor');
    %[t,s] = title(ax,['TCRs per Cluster: ',num2str(TPC)],['Mean MI: ', num2str(channel_sum),' bits']);
    %t.FontSize = 16;
    %s.FontAngle = 'italic';
    
    xticks(ax, xt);
    yticks(ax, yt);
    xtickformat(ax, xtf);
    ytickformat(ax, ytf);
    
    %if ~strcmp(mi_variable, "concentration")
     %   xtl = string(xticklabels(ax));
      %  xtl(2:2:end)=nan;
       % xticklabels(ax, xtl);
    %end

    %xlabel(ax, xlbl,'FontSize',16,'Interpreter','Latex')
    %ylabel(ax, ylbl,'FontSize',16,'Interpreter','latex')
    %ax.FontSize = 16;
    

    if savetrue
        f1 = gcf;
        fname = ['channel_cap_TPC',num2str(TPC),'.png'];
        saveas(f1, fname);
    end

end