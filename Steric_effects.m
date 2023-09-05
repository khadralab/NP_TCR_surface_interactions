%% Clustered Surface
clear all

% TCR Params (Clusters)
rSurf = 1000;
num_clusters = 15;
cluster_radius = 50;
tcr_per_cluster = 20;
num_tcr = 20*15;
tcr_radius = 5;

%{
% TCR Params (Uniform Surface)
rSurf = 1000;
num_clusters = 0;
cluster_radius = 50;
num_tcr = 20*15;
tcr_per_cluster = num_tcr / (rSurf^2 / cluster_radius^2 - num_clusters);
%}

tcr_params = [rSurf, num_clusters, cluster_radius, tcr_per_cluster, num_tcr];

% Kinetic Params
k0= 0.01;
kon=0.1;
koff=2;

kinetic_params = [k0, kon, koff];

% NP Params 
np_rho=1;
vh = 75;
np_radius = 100;

np_params = [vh, np_radius, np_rho];

% Cluster Surface
[nt_bound, nt_pd, aveBound, Bound_pdf] = NP_surface_binding(np_params, kinetic_params, tcr_params, 0);
aveCovered = [1:length(nt_bound)] * nt_bound' ./ sum(nt_bound);
nt = [1:length(nt_pd)];

figure()
subplot(211)
plot(nt, nt_pd, 'r-','Displayname','Pre-Bind')
hold on
plot(nt, nt_bound,'b-','DisplayName','Post-Bind')
xlabel('TCRs covered by NP')
ylabel('Probability')
xlim([1 50])
title(['Average Covered TCRs = ',num2str(aveCovered)]);
legend()

subplot(212)
plot(nt, Bound_pdf, 'r-','Displayname','Bound')
xlabel('TCRs Bound by NP')
ylabel('Probability')
title(['Average Bound TCRs = ',num2str(aveBound)]);
legend()

%--------------------------------------------------------------------------
% Steric hindrance
%--------------------------------------------------------------------------
disp('Starting here...........................................................')
Kav = [];
for i = nt
    K = np_avidity(vh,i,k0,kon,koff);
    Kav = [Kav,K];
end

% Probability of binding to surface as function of avidity (i.e. covered TCRs)
bind_dist = np_rho*Kav ./ (1 + np_rho * Kav);

figure()
plot(nt,bind_dist)
xlabel('Cov TCR')
ylabel('Bind Prob')

figure()
i=1;
epsilon = 0;
num_tcr = [num_tcr];
p_bind = [sum(nt_bound)];

epsilon = aveCovered / num_tcr(i);

aB = [aveBound];
aC = [aveCovered];

while epsilon < 1
    
    disp(sum(nt_pd))
    
    if mod(i,5) == 0
        
        aveCovered = [1:length(nt_pd)] * nt_pd' ./ sum(nt_pd);
        
        subplot(411)
        hold on
        plot(nt, nt_pd,'Displayname',[num2str(i),' | ',num2str(aveCovered)])
        %disp('Covered Pre')
        %disp(sum(nt_pd))

        aveCovered = [1:length(nt_bound)] * nt_bound' ./ sum(nt_bound);
        
        subplot(412)
        hold on
        plot(nt, nt_bound,'Displayname',[num2str(i),' | ', num2str(aveCovered)])
        %disp('Covered Post')
        %disp(sum(nt_bound))

        subplot(413)
        hold on
        plot(nt, p_bind(end)*Bound_pdf,'Displayname',num2str(i))
        %disp('Bound')
        %disp(sum(Bound_pdf))
    end
    
    nt_pd = nt_pd - epsilon * nt_bound;
    num_tcr = [num_tcr, num_tcr(i) - aveCovered];
    
    nt_bound = nt_pd .* bind_dist;
    p_bind = [p_bind, sum(nt_bound)];
    aveCovered = [1:length(nt_bound)] * nt_bound' ./ sum(nt_bound);
    epsilon = aveCovered / num_tcr(i);
    [aveBound, Bound_pdf] = boundTCRs(vh,max(nt),k0,kon,koff,nt_pd);
    
    aC = [aC, aveCovered];
    aB = [aB, aveBound];
    
    disp('Epsilon:');
    disp(epsilon);
    
    i=i+1;
end

subplot(411)
ylabel('Pdf')
title('Available TCRs (Pre-Bind)')
legend()

subplot(412)
ylabel('Pdf')
title('Covered TCRs (Post-Bind)')
legend()

subplot(413)
ylabel('Pdf')
xlabel('Number of TCRs')
title('Bound TCRs')
legend()

subplot(414)
%plot(1:i,num_tcr,'x--')
hold on
plot([1:i-1]-1, p_bind(1:end-1), 'x--')
ylabel('Probability of binding')
xlabel('Bound NPs')

figure()
subplot(221)
plot([0:i-2], aB(1:end-1))
xlabel('Bound NPs')
ylabel('Bound TCRs per NP')

subplot(222)
plot([0:i-2], aC(1:end-1))
xlabel('Bound NPs')
ylabel('Covered TCRs per NP')

subplot(223)
plot(cumsum(aB(1:end-1))/num_tcr(1), 1 - cumsum(aC(1:end-1))/num_tcr(1))
xlabel('Fraction of Bound Receptors')
ylabel('Fraction of Available Receptors')

%--------------------------------------------------------------------------
%% Local Functions
%
%
%% Compute the avidity of a given NP to the T cell surface
%--------------------------------------------------------------------------
function [Kav, Bound_pdf, Bound_tcr] = np_avidity(vh,nt,k0,kon,koff)
    if nt == 0
        Kav = 0;
        Bound_cdf = 0;
        Bound_tcr = 0;
        return
    end
    
    N = min(vh,nt);
    i = [1:N-1];
    f0 = k0*vh*nt;
    f = kon.*(vh-i).*(nt-i);
    f = [f0,f];
    b = [1:N] .* koff;
    
    
    y = cumprod(f./b);
    % Overall NP Avidity
    Kav = sum(y);
    
    % Distribution of Bound TCRs
    Y0 = 1 / (1 + Kav);
    Bound_pdf = [Y0, y.*Y0];
    Bound_cdf = cumsum(Bound_pdf);
    
    % Random number of bound TCRs estimated from above distribution
    r = rand(1);
    Bound_tcr = find(Bound_cdf >= r,1,'first')-1;
end

% Estimate bound TCRs
function [aveBound, bTCR_pd] = boundTCRs(vh,nt,k0,kon,koff,nt_bound)
    bTCR_pd=[];
    
    for i = 1:nt                  % Bound TCR

        if i > vh
            z =zeros(max(nt)-i+1,1)';
            bTCR_pd = [bTCR_pd, z];
            break
        end
        p_prod = 0;
        for j = i:nt            % Covered TCR
            p_cov = nt_bound(j) / sum(nt_bound);
            [Kav, Bound_pdf, Bound_tcr] = np_avidity(vh,j,k0,kon,koff);
            p_condi = Bound_pdf(i+1);
            p_prod = p_prod + p_condi * p_cov;
        end

        bTCR_pd = [bTCR_pd, p_prod];
    end
    bTCR_pd = bTCR_pd ./ sum(bTCR_pd);
    aveBound = [1:length(bTCR_pd)] * bTCR_pd';

end