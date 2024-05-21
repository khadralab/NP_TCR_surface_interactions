function y = PlotSurface(rSurf, tcr_pos, rTCR, np_pos, rNP)
    theta_vals = linspace(0,2*pi,1e3);
    
    
    f = figure('Color','white');
    plot(rSurf*cos(theta_vals), rSurf*sin(theta_vals),'-k')
    
    for i=1:size(tcr_pos,2)
        hold on
        plot(tcr_pos(1,i) + rTCR*cos(theta_vals), tcr_pos(2,i) + rTCR*sin(theta_vals),'-b')
    end
    
    for i=1:size(np_pos,2)
        hold on
        plot(np_pos(1,i) + rNP*cos(theta_vals), np_pos(2,i) + rNP*sin(theta_vals),'-r')
    end
    
    xlim([-rSurf-rNP rSurf+rNP])
    ylim([-rSurf-rNP rSurf+rNP])
    
    set(findall(f,'-property','FontSize'),'FontSize',16)
    set(findall(f,'-property','Interpreter'),'Interpreter','Latex')
end