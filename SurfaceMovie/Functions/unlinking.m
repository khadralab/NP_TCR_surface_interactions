function [tcr_states, bt, nt] = unlinking(tcr_states,b,e3,koff)
    cumb = cumsum(b) ./ sum(b);
    r1 = rand(1,e3);
    
    for i = 1:e3
        np_indx = find(cumb < r1(i), 1, 'last')+1;
        if isempty(np_indx)
            np_indx = 1;
        end
        tcr_indx = datasample(find(tcr_states(:,np_indx)>0),1,'Replace',false);
        tcr_states(tcr_indx, np_indx) = -1;
        
        bt = sum(tcr_states > 0, 1);
        nt = sum(tcr_states ~= 0, 1);
        b = koff*bt;
        cumb = cumsum(b) ./ sum(b);
        
        if sum(b) == 0
            return
        end
    end      
end