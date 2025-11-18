% Carrying capacity Monte Carlo sims
clear all

% TCR per Nanocluster:      0   |   3   |   5   |   10  |   20  |
% Nanocluster radii:        5   |50/s(6)|50/s(4)|50/s(2)|   50  |
% Number of clusters:       300 |   100 |   60  |   30  |   15  |
%
%


tpc = [1, 3, 5, 10, 20];
num_nc = [300, 100, 60, 30, 15];
tpN = [1, 3, 5, 5, 5];

nc_radius = [5, 50/sqrt(6), 50/sqrt(4), 50/sqrt(2), 50];
np_radius = [14,20];
tcr_radius = 5;

nsims = 1000;
attempts = 500;

tcr_per_np = zeros(1,nsims);

%{
for i = 1:nsims
    tcr_per_np(i) = MC_capacity(np_radius, tcr_radius, attempts);
end

mean_tpN = mean(tcr_per_np);
median_tpN = median(tcr_per_np);
max_tpN = max(tcr_per_np);
%}

mean_cap = zeros(length(np_radius), length(nc_radius));
med_cap = zeros(length(np_radius), length(nc_radius));
max_cap = zeros(length(np_radius), length(nc_radius));

for k = 1:length(np_radius)
    radius = np_radius(k) ./ nc_radius;
    
    capacity = zeros(length(radius), nsims);
    for i = 1:length(radius)
        disp('Loading...');
        disp([num2str(i / length(radius) * 100),' %']);
        r1 = 10;
        r2 = r1 * radius(i);

        for j=1:nsims
            cap = MC_capacity(r1,r2,attempts);     % r1 (domain) > r2 (steric exclusion)
            capacity(i,j) = cap; 
        end
    end

    mean_cap(k,:) = mean(capacity,2)';
    med_cap(k,:) = median(capacity,2)';
    max_cap(k,:) = max(capacity,[],2)';
    mode_cap(k,:) = mode(capacity,2)';

end

%%
f = figure('Color','White');
ax(1) = subplot(221);
ax(2) = subplot(222);
ax(3) = subplot(223);

hold(ax(:),'on')

col = parula(4);
for k = 1:length(np_radius)
    plot(ax(1),tpc, tpc ./ med_cap(k,:),':sq','Color', col(k,:),'DisplayName',[num2str(np_radius(k))],'LineWidth',1.5)
    plot(ax(2),tpc, med_cap(k,:),':sq','Color',col(k,:),'DisplayName',[num2str(np_radius(k))],'LineWidth',1.5);
    plot(ax(3),tpc, num_nc .* med_cap(k,:),':sq','Color',col(k,:),'DisplayName',[num2str(np_radius(k))],'LineWidth',1.5)
end
%plot(tpc, tpc' ./ med_cap, '-.xg','DisplayName','Median')
%plot(tpc, tpc' ./ max_cap, '--dr','DisplayName','Max')
l = legend(ax(1), 'Location','northeast');
xlabel(ax(:), 'TCRs per Cluster')
ylabel(ax(1),'TCRs per NP')
ylabel(ax(2),'Nanocluster NP Capacity')
ylabel(ax(3), 'Surface NP Capacity')

grid(ax(:),'on');
xticks(ax(:), tpc);
%yticks(gca,[1,2,3,4]);
title(l, 'NP Radius');

set(findall(f,'-property','FontSize'),'FontSize',16)
set(findall(f,'-property','Interpreter'),'Interpreter','Latex')
%% Local Functions

function [capacity] = MC_capacity(r1,r2, attempts)
    
    plt = false;
    pos=[];
    i=0;

    while i < attempts
        rho = rand(1,1) + rand(1,1);
        rho(rho>1)=2-rho(rho>1);
        theta = 2*pi*rand(1,1);

        % Candidate circle position
        temp_pos = [r1*rho.*cos(theta); r1*rho.*sin(theta)];

        % Check overlap
        d = dist(pos', temp_pos);
        mind = min(d);

        if mind < 2*r2
            i=i+1;
            continue
        else
            i=0;
            pos = [pos, temp_pos];
        end
    end

    capacity = size(pos,2);
    
    if plt
        figure('Color', 'White')
        angles = linspace(0,2*pi,1e3);

        plot(r1*cos(angles), r1*sin(angles), '-k');

        hold on

        for j = 1:size(pos,2)
            plot(r2*cos(angles)+pos(1,j), r2*sin(angles)+pos(2,j), '-b');
        end
    end
end