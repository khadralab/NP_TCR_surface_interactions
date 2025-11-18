clear all
addpath('Plotting/Functions/');
addpath('Functions/');

%% 20 TCR Cluster Surface
rSurf = 1000;
tcr_per_cluster = 20;
num_tcr = 300;
num_clusters = floor(num_tcr/tcr_per_cluster);
cluster_radius = 50;
rTCR = 5;
rNP = 20;

col = [1,0,0];

tcr_pos = gen_clusters(rSurf, cluster_radius, num_clusters, tcr_per_cluster, num_tcr);

PlotSurface(rSurf, tcr_pos, rTCR, [], rNP, col)

%% 10 TCR Cluster Surface
rSurf = 1000;
tcr_per_cluster = 10;
num_tcr = 300;
num_clusters = floor(num_tcr/tcr_per_cluster);
cluster_radius = 50/sqrt(2);
rTCR = 5;
rNP = 20;

tcr_pos = gen_clusters(rSurf, cluster_radius, num_clusters, tcr_per_cluster, num_tcr);

PlotSurface(rSurf, tcr_pos, rTCR, [], rNP, tcr_per_cluster)

%% 5 TCR Cluster Surface
rSurf = 1000;
tcr_per_cluster = 5;
num_tcr = 300;
num_clusters = floor(num_tcr/tcr_per_cluster);
cluster_radius = 50/2;
rTCR = 5;
rNP = 20;

tcr_pos = gen_clusters(rSurf, cluster_radius, num_clusters, tcr_per_cluster, num_tcr);

PlotSurface(rSurf, tcr_pos, rTCR, [], rNP, tcr_per_cluster)

%% 3 TCR Cluster Surface
rSurf = 1000;
tcr_per_cluster = 3;
num_tcr = 300;
num_clusters = floor(num_tcr/tcr_per_cluster);
cluster_radius = 50/sqrt(6);
rTCR = 5;
rNP = 20;

col = [0,0,1];

tcr_pos = gen_clusters(rSurf, cluster_radius, num_clusters, tcr_per_cluster, num_tcr);

PlotSurface(rSurf, tcr_pos, rTCR, [], rNP, tcr_per_cluster, col)
%% Uniform Surface
rSurf = 1000;
tcr_per_cluster = 0;
num_tcr = 300;
num_clusters = 20;
cluster_radius = 50;
rTCR = 5;
rNP = 20;

col = [0,0,0];

tcr_pos = gen_clusters(rSurf, cluster_radius, num_clusters, tcr_per_cluster, num_tcr);

PlotSurface(rSurf, tcr_pos, rTCR, [], rNP, col)
%% 20 TPC Mixed Surface

rSurf = 1000;
tcr_per_cluster = 20;
num_tcr = 300;
max_clusters = floor(num_tcr/tcr_per_cluster);
num_clusters = 10;
cluster_radius = 50;
rTCR = 5;
rNP = 20;

tcr_pos = gen_clusters(rSurf, cluster_radius, num_clusters, tcr_per_cluster, num_tcr);

col = [1, 0.6, 0];

PlotSurface(rSurf, tcr_pos, rTCR, [], rNP, col)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 20 TCR Single Cluster (Centered)

rSurf = 1000;
num_tcr = 20;
cluster_radius = 50;
rTCR = 5;
rNP = 20;

tcr_pos = generate_pos(cluster_radius, rTCR, num_tcr);
np_pos = generate_NPs(rSurf, 2*rNP, 1000, []);
[nt, d] = count_covered_tcrs(tcr_pos, np_pos, rNP);
np_pos = np_pos(:,nt>0);
nt = nt(nt>0);

PlotSurface(rSurf, tcr_pos, rTCR, np_pos, rNP)

%% 80 TCR Single Cluster (Centered)

rSurf = 1000;
cluster_radius = 100;
num_tcr = 20 * (cluster_radius / 50)^2;

rTCR = 5;
rNP = 20;

tcr_pos = generate_pos(cluster_radius, rTCR, num_tcr);

PlotSurface(rSurf, tcr_pos, rTCR, [], rNP)

%% 2000 TCR Single Cluster (Centered)

rSurf = 1000;
cluster_radius = 500;
num_tcr = 20 * (cluster_radius / 50)^2;

rTCR = 5;
rNP = 20;

tcr_pos = generate_pos(cluster_radius, rTCR, num_tcr);

PlotSurface(rSurf, tcr_pos, rTCR, [], rNP)
%% 8000 TCR Uniform Surface

rSurf = 1000;
tcr_per_cluster = 1;
num_tcr = 8000;
num_clusters = 0;
cluster_radius = 50;
rTCR = 5;
rNP = 20;

tcr_pos = gen_clusters(rSurf, cluster_radius, num_clusters, tcr_per_cluster, num_tcr);

PlotSurface(rSurf, tcr_pos, rTCR, [], rNP)