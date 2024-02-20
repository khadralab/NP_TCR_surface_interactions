function [tcr_states, np_pos] = updateBoundNPs(tcr_states, tcr_pos, np_pos, e1, rSurf, rNP, k0, vh)
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
    overlaps = dist_between_nps < 2*rNP;
    [I, ~] = find(overlaps);
    temp_np_pos(:,I) = [];
    
    [~, covMatrix] = count_covered_tcrs(tcr_pos, temp_np_pos, rNP);
    covMatrix = double(covMatrix);
    
    for i=1:size(temp_np_pos,2)
        K = find(np_pos(1,:) == 5000, 1);
        np_pos(:,K) = temp_np_pos(:,i);
        indx = datasample(find(covMatrix(i,:)==1),1,'Replace',false);
        covMatrix(i,indx)=-1;
        tcr_states(:,K) = -covMatrix(i,:)';
    end
    
    checkerror = sum(tcr_states~=0,2);

    if any(checkerror > 1)
        error('TCR duplications')
    end
end