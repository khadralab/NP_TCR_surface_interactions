function y = PlotSurface(ax, rSurf, tcr_pos, rTCR, np_pos, rNP,col)
    theta_vals = linspace(0,2*pi,1e3);
    
    
    %f = figure('Color','white');
    plot(ax, rSurf*cos(theta_vals), rSurf*sin(theta_vals),'-k','Linewidth',1.2)
    
    for i=1:size(tcr_pos,2)
        hold(ax,'on');
        plot(ax, tcr_pos(1,i) + rTCR*cos(theta_vals), tcr_pos(2,i) + rTCR*sin(theta_vals),'LineStyle','-','Color',col)
    end
    
    for i=1:size(np_pos,2)
        hold(ax,'on');
        plot(ax, np_pos(1,i) + rNP*cos(theta_vals), np_pos(2,i) + rNP*sin(theta_vals),'-r')
    end
    
    xlim(ax, [-rSurf-rNP rSurf+rNP])
    ylim(ax, [-rSurf-rNP rSurf+rNP])
    
    set(ax,'XTickLabel',[],'YTickLabel',[],'tickdir','out','box','off')
    axis(ax,'off');
    
    %set(findall(ax,'-property','FontSize'),'FontSize',16)
    %set(findall(ax,'-property','Interpreter'),'Interpreter','Latex')
    
    %fname = ['../Figures/Manuscript/Final-Figures/Surface-',num2str(tpc),'TPC'];
    %saveas(f, fname,'svg');
end