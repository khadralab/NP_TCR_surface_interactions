function [tcr_states] = crosslinking(tcr_states,f,e2,vh,kon)
    cumf = cumsum(f) ./ sum(f);
    r1 = rand(1,e2);
    
    if sum(tcr_states == -1,'all') < e2
        e2 = sum(tcr_states == -1,'all');
    end
    
    for i = 1:e2
        np_indx = find(cumf < r1(i), 1, 'last')+1;
        if isempty(np_indx)
            np_indx = 1;
        end
        
        if sum(f)==0
            return
        end
        
        if isempty(find(tcr_states(:,np_indx)==-1))
            error('Its empty!')
        end
        
        tcr_indx = datasample(find(tcr_states(:,np_indx)==-1),1,'Replace',false);
        tcr_states(tcr_indx, np_indx) = 1;
        
        bt = sum(tcr_states > 0, 1);
        nt = sum(tcr_states ~= 0, 1);
        f = kon*(vh-bt).*(nt-bt);
        cumf = cumsum(f) ./ sum(f);
    end  
end