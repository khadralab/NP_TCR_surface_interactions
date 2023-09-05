clear all
path(path, '/home/louis/Desktop/NP-Revival/Functions/')

%C = 1; U = 2;
choose_np = 1;
choose_surface = 1;

% Kinetic and NP Params

if choose_np == 1  
    k0= 0.1;
    kon=0.1;
    koff=0.1;
    
    np_rho=1;
    vh = 5;
    np_radius = 20;
    
elseif choose_np == 2
    k0= 0.01;
    kon=0.1;
    koff=1.22;
    
    np_rho=1;
    vh = 34;
    np_radius = 115;
    
else
    k0= 0.01;
    kon=0.1;
    koff=4.94;
    
    np_rho=0.02017;
    vh = 28.5;
    np_radius = 127;
end

% TCR Params
if choose_surface == 1
    rSurf = 3000;
    num_clusters = 15;
    cluster_radius = 50;
    tcr_per_cluster = 20;
    num_tcr = 20*15;
    tcr_radius = 5;
    
elseif choose_surface == 2
    rSurf = 3000;
    num_clusters = 0;
    cluster_radius = 50;
    num_tcr = 20*15;
    tcr_per_cluster = 20; %num_tcr / (rSurf^2 / cluster_radius^2 - num_clusters);
    tcr_radius = 5;
    
else
    rSurf = 1000;
    num_clusters = 15;
    cluster_radius = 50;
    tcr_per_cluster = 20;
    num_tcr = 20*15;
    tcr_radius = 5;
end

k0 = k0 * np_rho;

kinetic_params = [k0, kon, koff];
np_params = [vh, np_radius, np_rho];
tcr_params = [rSurf, num_clusters, cluster_radius, tcr_per_cluster, num_tcr];

% Distribution of covered TCRs and bound TCRs
%[nt_bound, nt_pd, bd_pd, Kav, p_bind] = NP_surface_binding(np_params, kinetic_params, tcr_params, 0);

% Surface Coverage and Selectivity
Selectivity(np_params, kinetic_params, tcr_params, 1);

% MC Surface Coverage
[K_np, p_bound, nt, BoundTCRs] = MC_Surface(np_params, kinetic_params, tcr_params, 1);

%%

vh = 1;
nt = 10;
k0= 0.01;
kon=0.1;

figure()

for nt = [1,5,10,50]
    
    kav = [];
    theta = [];

    for koff = [-1:0.1:1]
        koff = 10^koff;

        k = np_avidity(vh,nt,k0,kon,koff);
        kav = [kav, k];
        theta = [theta, k / (1+k)];
    end
    
    subplot(211)
    plot(10.^[-1:0.1:1], kav, 'DisplayName', [num2str(nt)])
    hold on
    ylabel('Avidity')
    xlabel('Koff')
    legend()
    
    subplot(212)
    plot(10.^[-1:0.1:1], theta, 'DisplayName', [num2str(nt)])
    hold on
    ylabel('Surface Coverage')
    xlabel('Koff')
    legend()
end



%%


cluster_params = [3000, 15, 50, 20, 15*20];
uni_params = [3000,0,50,20,15*20];
kinetic_params = [0.01,0.1,20];

% NP Params (v, r, rho)
weak_np = [50, 100, 1];
strong_np = [50, 100, 1];

weak_select = [];
strong_select = [];
cluster_strong = [];
cluster_weak = [];
uni_weak = [];
uni_strong = [];

x = -4:0.5:0;
x = 10.^x;

for rho = x
    
    weak_kinetics = [rho*0.01,0.1,5];
    strong_kinetics = [rho*0.01, 0.1, 2];
    
    [selectivity, p_bound] = Selectivity(weak_np, weak_kinetics, cluster_params, 0);
    weak_select = [weak_select, find(selectivity == max(selectivity))];
    
    [selectivity, p_bound] = Selectivity(strong_np, strong_kinetics, cluster_params, 0);
    strong_select = [strong_select, find(selectivity == max(selectivity))];
    
    [nt_bound, nt_pd, aveBound, bd_pd] = NP_surface_binding(strong_np, strong_kinetics, cluster_params, 0);
    cluster_strong = [cluster_strong, aveBound];
    
    [nt_bound, nt_pd, aveBound, bd_pd] = NP_surface_binding(strong_np, strong_kinetics, uni_params, 0);
    uni_strong = [uni_strong, aveBound];
    
    [nt_bound, nt_pd, aveBound, bd_pd] = NP_surface_binding(weak_np, weak_kinetics, cluster_params, 0);
    cluster_weak = [cluster_weak, aveBound];
    
    [nt_bound, nt_pd, aveBound, bd_pd] = NP_surface_binding(weak_np, weak_kinetics, uni_params, 0);
    uni_weak = [uni_weak, aveBound];
end

figure()
plot(x, weak_select,'DisplayName','Weak NP')
hold on
plot(x, strong_select,'DisplayName','Strong NP')
xlabel('NP Concentration')
ylabel('Number of TCRs of peak selectivity')
legend()

figure()
semilogy(x, cluster_strong, 'r-', 'DisplayName','Clus-Strong')
hold on
semilogy(x, cluster_weak, 'r--', 'DisplayName','Clus-Weak')
semilogy(x, uni_strong, 'b-','DisplayName', 'Uni-Strong')
semilogy(x, uni_weak,'b--','DisplayName','Uni-Weak')
xlabel('NP Concentration')
ylabel('Bound TCRs')
legend()

%% 
rSurf = 3000;
num_clusters = 15;
cluster_radius = 50;
tcr_per_cluster = 20;
num_tcr = 20*15;
tcr_radius = 5;

np_rho=1;
vh = 5;
np_radius = 20;

np_params = [vh, np_radius, np_rho];
tcr_params = [rSurf, num_clusters, cluster_radius, tcr_per_cluster, num_tcr];

i=1;
cnp=zeros(500,1);
while i < 500
    cnp(i) = CarryingCapacity(np_params, tcr_params,0);
    i=i+1;
end
mean_cnp = mean(cnp);

figure()
histogram(cnp,'Normalization','pdf')
xlabel('Maximum Bound NPs')
ylabel('PDF')
title('Carrying Capacity')

fname = ['Cluster_cc'];
saveas(gca,fname, 'png');
%% Compute the avidity of a given NP to the T cell surface
%--------------------------------------------------------------------------
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