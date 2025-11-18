function mi_cell = Valence_MI_matrix(ifsave, norm)

    % New Simulations
    koff_vals = [0.0001, 0.001, 0.01, 0.1];
    
    % NP Params
    np_radius = 20;
    %np_valence = 5;
    valence_vals = [1,2,3,4,5,10];
    rho_vals = linspace(-4,3.5,16);
    rho_vals = 10.^ rho_vals;
    
    % Surface Params
    surfaces = [1,3,5,10,20];
    
    mi_cell = cell(2, length(surfaces));
    channel_cap=ones(1,length(surfaces));
    
    for t = 1:length(surfaces)
        TPC = surfaces(t);
    
        mi = zeros(length(rho_vals), length(koff_vals));
        variance_mi = zeros(length(rho_vals), length(koff_vals));
    
        for i = 1:length(rho_vals)
            np_rho = rho_vals(i);
            for j = 1:length(koff_vals)
                koff = koff_vals(j);
                tcr_pdf = [];
                k=1;
                for np_valence = valence_vals
                    [bound_tcr, bound_np] = load_BoundTCR_dist(koff, np_rho, np_radius, np_valence, TPC);
                    id = k*ones(size(bound_tcr));
                    tcr_pdf = [tcr_pdf; [bound_tcr, id]];
    
                    k=k+1;
                end
            
                tcr_dist = tcr_pdf(:,1);
                v_dist = tcr_pdf(:,2);
            
                [mean_mi, var_mi] = bootstrap_MI(tcr_dist, v_dist, 30, norm);
                mi(i,j) = mean_mi;
                variance_mi(i,j) = var_mi;
            end
        end
        
        mi_cell{1,t} = mi;
        mi_cell{2,t} = variance_mi;
        channel_cap(t) = max(mi,[],'all');
    
    end

    if ifsave && norm
        save('Mutual Information/valence_nmi','mi_cell');
    elseif ifsave && ~norm
        save('Mutual Information/valence_mi', 'mi_cell');
    end
end