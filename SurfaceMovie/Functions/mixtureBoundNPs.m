function [tcr_states, np_pos, type_array] = mixtureBoundNPs(tcr_states, tcr_pos, np_pos, np_type, type_array, e1, rSurf, np1_params, np2_params)
    % NP 1 params
    rNP1 = np1_params(1);
    vh1 = np1_params(2);
    k01 = np1_params(3);
    
    % NP 2 params
    rNP2 = np2_params(1);
    vh2 = np2_params(2);
    k02 = np2_params(3);
    
    if np_type == 1
        rNP = rNP1;
        vh = vh1;
        k0 = k01;
    else
        rNP = rNP2;
        vh = vh2;
        k0 = k02;
    end


    free_tcr_pos = tcr_pos(:,sum(tcr_states,2) == 0);
    
    checkerror = sum(tcr_states~=0,2);

    if any(checkerror > 1)
        error('TCR duplications')
    end
    
    if isempty(free_tcr_pos)
        return
    end
    
    if e1 > 200
        e1 = 200;
    end
    
    temp_np_pos = generate_NPs(rSurf, 2*rNP, e1, []);
    
    [temp_nt, ~] = count_covered_tcrs(free_tcr_pos, temp_np_pos, rNP);
    %temp_nt = sum(tcr_indx);
    temp_np_pos = temp_np_pos(:,temp_nt > 0);
    temp_nt = temp_nt(temp_nt > 0);
    
    if isempty(temp_np_pos)
        return
    end
    
    % Probability of new NP binding
    pbind = k0*vh*temp_nt ./ (1 + k0*vh*temp_nt);
    r1 = rand(1,length(pbind));
    
    % Filter successful bindings
    temp_np_pos = temp_np_pos(:, r1 < pbind);
        
    if isempty(temp_np_pos)
        return
    end
    
    % Check for overlap with other NPs
    dist_between_nps = dist(temp_np_pos', np_pos);
    radius_array = rNP * ones(size(dist_between_nps));
    cols = rNP1 * (type_array==1) + rNP2 *(type_array==2);
    radius_array = radius_array + cols;
    overlaps = dist_between_nps < radius_array;
    [I, ~] = find(overlaps);
    temp_np_pos(:,I) = [];
    
    [~, covMatrix] = count_covered_tcrs(tcr_pos, temp_np_pos, rNP);
    covMatrix = double(covMatrix);
    
    for i=1:size(temp_np_pos,2)
        K = find(np_pos(1,:) == 5000, 1);
        np_pos(:,K) = temp_np_pos(:,i);
        type_array(K) = np_type;
        indx = datasample(find(covMatrix(i,:)==1),1,'Replace',false);
        covMatrix(i,indx)=-1;
        tcr_states(:,K) = -covMatrix(i,:)';
    end
    
    checkerror = sum(tcr_states~=0,2);

    if any(checkerror > 1)
        error('TCR duplications')
    end
end