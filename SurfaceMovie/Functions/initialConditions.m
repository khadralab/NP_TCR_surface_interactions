function [tcr_states, np_pos, bound_tcr, covered_tcr] = initialConditions(tcr_states, tcr_pos, np_pos, rSurf, rNP, vh)
    maxNP = 100;
    countNPs = floor(maxNP*rand(1));
    boundNPs = 0;
    
    iter=0;
    while (boundNPs < countNPs) && (iter < 500)
        [tcr_states, np_pos] = updateBoundNPs(tcr_states, tcr_pos, np_pos, countNPs, rSurf, rNP, 1000, vh);
        
        boundNPs = sum(np_pos(1,:) ~= 5000);
        iter = iter + 1;
    end
    
    tcr_states = arrayfun(@randstate, tcr_states);
    
    bound_tcr = sum(tcr_states > 0, 1);         % Bound TCRs per NP
    covered_tcr = sum(tcr_states ~= 0, 1);      % Covered TCRs per NP 
end

function [y] = randstate(x)
    if x == 0
        y = 0;
        return
    elseif x == -1
        y = randsample([-1,1:3],1);
        return
    else
        y = randsample([1:3],1);
        return
    end
end