%% Comprehensive Dynamics of the Contact Area Model
%
%


%% NP Avidity: Dependence on vh, nt and koff

k0 = 0.1;
kon = 2.5;
koff = 0.1;

figure()

for vh = [1,2,3,4,5]

    avidity = [];

    for nt = 1:100
        Kav = np_avidity(vh,nt,k0,kon,koff);
        avidity = [avidity, Kav];
    end

    subplot(211)
    semilogy([1:100], avidity, 'DisplayName',['vh = ',num2str(vh)],'Linewidth',1.5)
    hold on
    
    nt = 5;
    avidity = [];
    for koff = 1:100
        Kav = np_avidity(vh,nt,k0,kon,koff);
        avidity = [avidity, Kav];
    end
    
    subplot(212)
    loglog([1:100], avidity, 'DisplayName',['vh = ',num2str(vh)],'Linewidth',1.5)
    hold on
end
subplot(211)
xlabel('Covered TCRs')
ylabel('NP Avidity')
legend('Location','southeast','Fontsize',10)
set(gca,'Fontsize',12)
grid on

subplot(212)
xlabel('Koff')
ylabel('NP Avidity')
set(gca,'Fontsize',12)
grid on

%% Range of TCR densities

valences = [50,50];                                                   % Range of NP valences
dissociation = [1,100];                                  % Range of f0 affinities to ensure identical half-max activation for each valence
tcr_density = [1:50];
K_av = zeros(length(valences),length(tcr_density));

for i = 1:length(dissociation)
    vh = valences(i);
    k0 = 0.01;
    kon = 0.1;
    koff= dissociation(i);
    for nt = tcr_density
        N = min(vh,nt);
    
        j = [1:N-1];
        f0 = k0*vh*nt;
        f = kon.*(vh-j).*(nt-j);
        f = [f0,f];
        b = [1:N] .* koff;

        avidity = cumprod(f./b);
        K_av(i,nt) = sum(avidity);
    end
end
bound_fraction = K_av ./ (1+K_av);

labels = ['v = '] + string(valences)+['; Koff = ']+ string(dissociation);
%labels = {['Monovalent Strong'], ['Monovalent Weak'],['Multivalent Weak']};

figure()
plot(tcr_density,bound_fraction,'Linewidth',2);
grid on
leg = legend(labels,'Location','northwest');
xlabel('Number of Receptors')
ylabel('Fraction of bound sites')
%title(leg, 'Valence-log(Ka)')
set(gca,'Fontsize',16)

figure()
semilogx(tcr_density,bound_fraction.*tcr_density);
leg = legend(labels,'Location','southeast');
xlabel('Number of Receptors')
ylabel('Number of bound receptors')
title(leg, 'Valence-log(Ka)')

%% TCR Landscapes and Distributions of Covered TCRs
% Generate different types of clusters and TCR landscapes. Compute
% distributions of covered TCRs.
clear all 

% Kinetic
kon = 0.1;
k0 = 0.01;
koff= 20;

% NP
vh = 15;
rNP = 100;
np_rho = 1;
np_num = 1000;

% TCR NC
rTCR = 5;
cluster_radius = 50;
tcr_per_cluster = 20;
rSurf = 1000;
max_clusters = 15;
num_tcr = max_clusters * tcr_per_cluster;

for num_clusters = [15,0]
    disp(['Clusters = ',num2str(num_clusters)])
    
    tcr_pos = gen_clusters(rSurf, cluster_radius, num_clusters, tcr_per_cluster, num_tcr);    
    %tcr_pos = generate_pos(rSurf,rTCR,50*tcr_num);
    %tcr_pos = generate_clusters(tcr_pos, rSurf, num_clusters, cluster_radius,tcr_num, epsilon);
    
    num_NP = 100000;
    rho = rand(1,num_NP) + rand(1,num_NP);
    rho(rho>1)=2-rho(rho>1);
    theta = 2*pi*rand(1,num_NP);

    np_pos = [rSurf*rho.*cos(theta); rSurf*rho.*sin(theta)];
    
    nt = covered_tcrs(tcr_pos, np_pos, rNP);
    
    K_np = [];
    p_bound = [];
    for i=1:size(np_pos,2)
        avidity = np_avidity(vh,nt(i),k0,kon,koff);
        K_np = [K_np, avidity];
        p_bound = [p_bound, np_rho.*avidity /(1 + np_rho.*avidity)];
    end
    
    disp('Avidity determined...')
    
    x = rand(1,length(np_pos));
    np_pos = np_pos(:,x<p_bound);
    nt_bound = nt(x<p_bound);
    angles = linspace(0,2*pi,500);
    
    
    figure()
    subplot(2,2,1:2)
    angles = linspace(0,2*pi,500);
    plot(rSurf*cos(angles), rSurf*sin(angles),'k-','Linewidth',2)
    hold on
    plot(rNP*cos(angles), rNP*sin(angles),'r-','Linewidth',1)
    for i = 1:length(tcr_pos)
        plot(rTCR*cos(angles)+tcr_pos(1,i),rTCR*sin(angles)+tcr_pos(2,i),'b-')
    end
    xlim([-rSurf rSurf]);
    ylim([-rSurf rSurf]);
    dim = [0.1 0.6 0.3 0.3];
    axis off
    
    % Plot histograms with fits of each distribution
    
    nt(nt==0) = [];
    pd1 = fitdist(nt','gamma');
    pd2 = fitdist(nt','poisson');
    bins = 0:max(nt);
    y1 = pdf(pd1,bins);
    y2 = pdf(pd2,bins);
    
    bin_lim = round(max(nt));
    
    subplot(2,2,3:4)
    histogram(nt,(0.5:bin_lim + 0.5), 'Normalization','probability','DisplayName','MC Sim')
    hold on
    %plot(bins,y1,'r-')
    %plot(bins,y2,'g-')
    xlabel('Covered TCRs (exc. nt=0)')
    ylabel('Frequency')
    %title('MC Simulations')
    
    plt = false;
    nt = coveredTCRs(rSurf, num_tcr, num_clusters, cluster_radius, tcr_per_cluster, rNP, plt);
    nt(nt==0)=[];
    pd = pdf_CovTCR(rSurf, num_tcr, num_clusters, cluster_radius, tcr_per_cluster, rNP);
    
    %subplot(224)
    %histogram(nt,(0.5:bin_lim + 0.5),'Normalization','probability')
    hold on
    plot(bins+1, pd(2:length(bins)+1)./sum(pd(2:length(bins)+1)),'r-', 'DisplayName','Theoretical','Linewidth',2)
    %plot(bins,y1,'r-')
    %plot(bins,y2,'g-')
    %xlabel('Covered TCRs (exc. nt=0)')
    %ylabel('Frequency')
    %title('Theoretical')
    legend()
    
end

%% Binding of NPs on different cell surfaces
% Compare the number of bound NPs for different T cell surface
% organizations.
plt = true;

% Kinetic
kon = 0.1;
k0 = 0.01;
koff= 53;

% NP
vh = 101.8;
rNP = 71.6;
np_rho = 0.14;
np_num = 1000;

% TCR NC
cluster_radius = 50;
tcr_per_cluster = 20;
max_clusters = 15;
num_tcr = max_clusters * tcr_per_cluster;

 % Surface geometry
rSurf = 1000;
rTCR = 5;

for num_clusters = [15,10,0]
    
    disp(['Clusters:', num2str(num_clusters)])
   
    tcr_pos = gen_clusters(rSurf, cluster_radius, num_clusters, tcr_per_cluster, num_tcr);  
    disp('TCR positions generated...')
    np_pos = generate_pos(rSurf, rNP, np_num);
    disp('NP positions generated...')
    nt = covered_tcrs(tcr_pos,np_pos,rNP);
    disp('Covered TCRs computed...')

    K_np = [];
    p_bound = [];
    for i=1:size(np_pos,2)
        avidity = np_avidity(vh,nt(i),k0,kon,koff);
        K_np = [K_np, avidity];
        p_bound = [p_bound, np_rho.*avidity /(1 + np_rho.*avidity)];
    end
    disp('Avidity determined...')
    x = rand(1,length(np_pos));
    np_pos = np_pos(:,x<p_bound);
    angles = linspace(0,2*pi,500);
    
    tot_covTCR = length(np_pos) * mean(nt(nt~=0));


    %DB_SCAN(tcr_pos',5)
    %disp('DBSCAN Completed...')

    % Plotting Monte Carlo Sims
    if plt == true
        figure()
        plot(rSurf*cos(angles), rSurf*sin(angles),'k-','Linewidth',2)
        hold on
        for i = 1:length(tcr_pos)
            plot(rTCR*cos(angles)+tcr_pos(1,i),rTCR*sin(angles)+tcr_pos(2,i),'b-')
        end
        for i = 1:size(np_pos,2)
            plot(rNP*cos(angles)+np_pos(1,i), rNP*sin(angles)+np_pos(2,i),'r-')
        end
        
        h = zeros(3, 1);
        h(1) = plot(nan,nan,'or');
        h(2) = plot(nan,nan,'ob');
        h(3) = plot(nan,nan,'ok');
        
        xlim([-rSurf rSurf]);
        ylim([-rSurf rSurf]);
        dim = [0.1 0.6 0.3 0.3];
        str = {['Bound NPs: ',num2str(length(np_pos))],['Average number of TCRs covered per NP: ',num2str(mean(nt(nt~=0)))],['Total Covered TCRs: ',num2str(tot_covTCR)]};
        %annotation('textbox',dim,'String',str,'FitBoxToText','on');
        legend(h,str,'location','southoutside','Fontsize', 12)
        axis off
    end
end

%% Distributions of covered TCRs by Bound NPs
clear all

% TCR Params
rSurf = 1000;
num_clusters = 15;
cluster_radius = 50;
tcr_per_cluster = 20;
num_tcr = 20*15;

% Kinetic Params
k0= 0.01;
kon=0.1;
koff=5;

% NP Params
np_rho=1;
vh = 150;
np_radius = 150;

nt = [1:50];
nt_pd = pdf_CovTCR(rSurf, num_tcr, num_clusters, cluster_radius, tcr_per_cluster, np_radius);
nt_pd(1) = [];
nt_pd(length(nt)+1:end)=[];
Kav = [];
for i = nt
    K = np_avidity(vh,i,k0,kon,koff);
    Kav = [Kav,K];
end
p_bind = np_rho*Kav ./ (1 + np_rho * Kav);

nt_bound = nt_pd .* p_bind;
bind_freq = sum(nt_bound);

figure()
subplot(211)
plot(nt, nt_pd, 'r-','Displayname','Pre-Bind')
hold on
plot(nt, nt_bound, 'g-','Displayname','Post-Bind')
xlabel('TCRs covered by NP')
ylabel('Probability')
legend()
title(['Cluster Landscape: P_{bind} =',num2str(bind_freq)])

subplot(212)
plot(nt, p_bind, 'b-','Displayname','Bind Prob.')
xlabel('TCRs covered by NP')
ylabel('Probability')
legend()

% TCR Params
rSurf = 1000;
num_clusters = 20*15;
cluster_radius = 50;
tcr_per_cluster = 1;
num_tcr = 20*15;

nt_pd = pdf_CovTCR(rSurf, num_tcr, num_clusters, cluster_radius, tcr_per_cluster, np_radius);
nt_pd(1) = [];
nt_pd(length(nt)+1:end)=[];

nt_bound = nt_pd .* p_bind;
bind_freq = sum(nt_bound);

figure()
subplot(211)
plot(nt, nt_pd, 'r-','Displayname','Pre-Bind')
hold on
plot(nt, nt_bound, 'g-','Displayname','Post-Bind')
xlabel('TCRs covered by NP')
ylabel('Probability')
legend()
title(['Uniform Landscape: P_{bind}:', num2str(bind_freq)])

subplot(212)
plot(nt, p_bind, 'b-','Displayname','Bind Prob.')
xlabel('TCRs covered by NP')
ylabel('Probability')
legend()

%%  Selectivity

clear all

% TCR Params
rSurf = 1000;
num_clusters = 15;
cluster_radius = 50;
tcr_per_cluster = 20;
num_tcr = 20*15;

% Kinetic Params
k0= 0.01;
kon=0.1;
koff=5;

% NP Params
np_rho=1;
vh = 10;
np_radius = 100;

nt = [1:150];

K_np = [];
p_bound = [];
selectivity = [];

for i=nt
    avidity = np_avidity(vh,i,k0,kon,koff);
    K_np = [K_np, avidity];
    p_bound = [p_bound, np_rho.*avidity /(1 + np_rho.*avidity)];
    
    if i>1
        alpha = log10(p_bound(i))-log10(p_bound(i-1));
        alpha = alpha / (log10(i)-log10(i-1));
        selectivity = [selectivity,alpha];
    end
end

x = nt ./ (pi * np_radius^2);
x1 = tcr_per_cluster ./ (pi * cluster_radius^2);
x2 = num_tcr ./ (pi * rSurf^2);

figure()
subplot(211)
plot(x, p_bound)
xline(x1,'g--')
xline(x2,'r--')
ylabel('Bind Prob')

subplot(212)
plot(x(1:end-1), selectivity)
xline(x1,'g--')
xline(x2,'r--')
xlabel('TCR Density')
ylabel('Selecivity')

%% Plot MC Results

figure()
cols = parula(np_num);

for i = 1:np_num
    subplot(311)
    plot(mean(nt,1))
    hold on

    subplot(312)
    plot(mean(Kav,1))
    hold on
    
    subplot(313)
    errorbar(mean(p_bind,1),var(p_bind,1))
    hold on
end

subplot(311)
xlabel('Covered TCRs')


subplot(312)
xlabel('NP Avidity')


subplot(313)
xlabel('Probability of binding')

%% Local Functions
%
%
%% Spherical contact distribution
function [cdf, pdf] = SCD(x, lambda)
    cdf = 1-exp(-lambda * pi .* (x.^2));
    pdf = 2*lambda*pi.*x.*exp(-lambda*pi.*(x.^2));
end

%% Area of Overlap
function area = overlap_area(x, rnc, rnp)
% x =  NP distance from center of NC
% rnc = Radius of nanocluster
% rnp = Radius of nanoparticle

    if rnp > rnc
        tmp = rnc;
        rnc = rnp;
        rnp = tmp;
    end

    h = 1./(2.*x).*sqrt(2.*x.^2.*rnp^2+2.*x.^2.*rnc^2+2*rnp^2*rnc^2-rnp^4-rnc^4-x.^4);
    
    a1 = rnp^2.*asin(h./rnp)+rnc^2.*asin(h./rnc)-x.*h;
    
    a2 = rnp^2 .*asin(h./rnp)-h.*sqrt(rnp^2-h.^2);
    a2 = a2 - rnc^2 .*asin(h./rnc)+ h.*sqrt(rnc^2-h.^2);
    a2 = pi * rnp^2 - a2;
    
    area = a1;
    
    cond = sqrt(rnc^2 -h.^2);
    
    area(x <= cond) = a2(x <= cond);
    
    area = real(area);

end
%% Poisson rate parameter for average covered TCRs
function lambda = ave_CovTCR(x, rnc, rnp, tcr_rho, free_rho)
    A = overlap_area(x, rnc, rnp);
    lambda = A .* tcr_rho + (pi * rnp^2 - A).*free_rho;
end

%% Generating random positions of points on disk
function pos = generate_pos(rSurf, dmin, num_points);
    rho = rand(1,num_points) + rand(1,num_points);
    rho(rho>1)=2-rho(rho>1);
    theta = 2*pi*rand(1,num_points);

    pos = [rSurf*rho.*cos(theta); rSurf*rho.*sin(theta)];

    i=2;

    while i<=length(pos)
        d = dist(pos(:,i)',pos(:,1:i-1));
        min_d = min(d);
        nsims = 0;
        while min_d < 2*dmin && nsims < 500
            rho = rand(1,1)+rand(1,1);
            rho(rho>1) = 2-rho(rho>1);
            th = 2*pi*rand(1,1);
            pos(:,i) = [rSurf*rho*cos(th); rSurf*rho*sin(th)];
        
            d = dist(pos(:,i)',pos(:,1:i-1));
            min_d = min(d);
            nsims = nsims+1;
        end
        
        if nsims == 500
            pos(:,i)=[];
            i=i-1;
        end
        i=i+1;
    end
end
%
%
%
%% Generate TCR Nanoclusters (Alternative Gaussian Filtering)
function tcr_pos = gen_clusters(rSurf, cluster_radius, num_clusters, tcr_per_cluster, num_tcr)
    rTCR = 5;
    nc_pos = generate_pos(rSurf-cluster_radius, 2*cluster_radius, num_clusters);
    tcr_pos = [];
    for i = 1:num_clusters
        %lambda = poissrnd(tcr_per_cluster);
        lambda = tcr_per_cluster;
        temp_pos = generate_pos(cluster_radius, rTCR, lambda);
        temp_pos = temp_pos + nc_pos(:,i);
        tcr_pos = [tcr_pos, temp_pos];
    end
    
    free_tcr = num_tcr -length(tcr_pos);
    
    if free_tcr < 0
        tcr_pos = tcr_pos(:,1:num_tcr);
    else
        temp_pos = generate_pos(rSurf, rTCR, free_tcr);
        tcr_pos = [tcr_pos, temp_pos];
    end
end
%% Compute the number of covered TCRs for each NP
function nt = covered_tcrs(tcr_pos, np_pos, rNP)
    nt = zeros(1,size(np_pos,2));
    for i = 1:size(np_pos,2)
        d = dist(np_pos(:,i)',tcr_pos);
        nt(i) = length(d(d<rNP));
    end
end
%
%
%
%% Compute the avidity of a given NP to the T cell surface
function Kav = np_avidity(vh,nt,k0,kon,koff)
    if nt == 0
        Kav = 0;
        return
    end
    
    N = min(vh,nt);
    i = [1:N-1];
    f0 = k0*vh*nt;
    f = kon.*(vh-i).*(nt-i);
    f = [f0,f];
    b = [1:N] .* koff;
    
    Kav = cumprod(f./b);
    Kav = sum(Kav);
end
%
%
%
%
%% Use DBSCAN to perform clustering on TCR positions and yield cluster labels.
function labels = DB_SCAN(pos, minpts)
    kD = pdist2(pos,pos,'euc','Smallest',minpts);
    
    figure()
    plot(sort(kD(end,:)));
    title('k-distance graph')
    xlabel(['Points sorted with ',num2str(minpts),'th nearest distances'])
    ylabel([num2str(minpts),'th nearest distances'])
    grid
    
    epsilon = 40;
    
    labels = dbscan(pos,epsilon,minpts);
    
    figure()
    numGroups = length(unique(labels));
    gscatter(pos(:,1),pos(:,2),labels,hsv(numGroups));
    title(['epsilon = ',num2str(epsilon),' and minpts = ',num2str(minpts)])
    grid
    
    sizeGroups = [];
    
    for  i = 1:numGroups
        sG = length(labels(labels == i));
        sizeGroups = [sizeGroups, sG];
    end
    
    figure()
    histogram(sizeGroups, 10)
    xlabel('TCRs per Nanocluster')
    ylabel('Count')
end

%% Cell Surface Binding
function [K_np, np_pos, nt] = surface_binding(np_params, tcr_params, kinetic_params, plt)
% NP Params = [ NP Radius, Effective Valence, Number of NP Positions, NP concentration]
% TCR Params = [ TCR Radius, TCR Density, Number of Clusters, Radius of Clusters]
% Kinetic Params = [ k0, kon, koff]
    
    rSurf = 1000;
    
    % Kinetic Parameters
    k0 = kinetic_params(1);
    kon = kinetic_params(2);
    koff = kinetic_params(3);
    
    % TCR Properties
    rTCR = tcr_params(1);
    tcr_density = tcr_params(2) / (50^2);
    tcr_num = round(tcr_density * rSurf^2);
    num_clusters = tcr_params(3);
    cluster_radius = tcr_params(4);
    epsilon = tcr_params(5);

    % NP Properties
    rNP = np_params(1);
    vh = np_params(2);
    np_num= np_params(3);
    np_rho= np_params(4);
    
    % Generate positions of TCRs, Clusters and NPs. Estimate covered TCRs.
    tcr_pos = generate_pos(rSurf,rTCR,20*tcr_num);
    tcr_pos = generate_clusters(tcr_pos, rSurf, num_clusters, cluster_radius,tcr_num, epsilon);
    disp('TCR positions generated...')
    np_pos = generate_pos(rSurf, rNP, np_num);
    disp('NP positions generated...')
    nt = covered_tcrs(tcr_pos,np_pos,rNP);
    disp('Covered TCRs computed...')

    % Compute avidity of each NP and probability of binding.
    K_np = [];
    p_bound = [];
    for i=1:size(np_pos,2)
        avidity = np_avidity(vh,nt(i),k0,kon,koff);
        K_np = [K_np, avidity];
        p_bound = [p_bound, np_rho.*avidity /(1 + np_rho.*avidity)];
    end
    disp('Avidity determined...')
    x = rand(1,length(np_pos));
    np_pos = np_pos(:,x<p_bound);
    nt = nt(x<p_bound);
    K_np = K_np(x<p_bound);
    angles = linspace(0,2*pi,500);


    %DB_SCAN(tcr_pos',5)
    %disp('DBSCAN Completed...')

    % Plotting Monte Carlo Sims
    if plt == true
        figure()
        plot(rSurf*cos(angles), rSurf*sin(angles),'k-','Linewidth',2)
        hold on
        for i = 1:length(tcr_pos)
            plot(rTCR*cos(angles)+tcr_pos(1,i),rTCR*sin(angles)+tcr_pos(2,i),'b-')
        end
        for i = 1:size(np_pos,2)
            plot(rNP*cos(angles)+np_pos(1,i), rNP*sin(angles)+np_pos(2,i),'r-')
        end
        
        h = zeros(3, 1);
        h(1) = plot(nan,nan,'or');
        h(2) = plot(nan,nan,'ob');
        h(3) = plot(nan,nan,'ok');
        
        xlim([-rSurf rSurf]);
        ylim([-rSurf rSurf]);
        dim = [0.1 0.6 0.3 0.3];
        str = {['Bound NPs: ',num2str(length(np_pos))],['Average number of TCRs covered per NP: ',num2str(mean(nt(nt~=0)))],['Effective valence of NPs: ',num2str(vh)]};
        %annotation('textbox',dim,'String',str,'FitBoxToText','on');
        legend(h,str,'location','southoutside','Fontsize', 12)
        axis off
    end
end