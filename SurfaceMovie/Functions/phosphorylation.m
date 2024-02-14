function [tcr_states] = phosphorylation(tcr_states, e4, maxP)
    [I, J] = find(tcr_states > 0 & tcr_states < maxP);
    
    if length(I) < e4
        e4 = length(I);
    end
    
    indx = datasample([1:length(I)],e4,'Replace',false);
    
    for i = 1:e4
        tcr_states(I(i), J(i)) = tcr_states(I(i), J(i)) + 1;
    end
end