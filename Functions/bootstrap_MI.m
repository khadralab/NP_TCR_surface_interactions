function [mean_mi, var_mi] = bootstrap_MI(X,Y,nstraps,norm)
    data = [X,Y];
    k = floor(0.2 * length(X));         % Sample 20% of the data.

    mi_array = zeros(1, nstraps);
    
    for i = 1:nstraps
        bootsample = datasample(data,k);

        bootX = bootsample(:,1);
        bootY = bootsample(:,2);

        mi_array(i) = mutual_info(bootX, bootY, norm);
    end

    mean_mi = mean(mi_array);
    var_mi = var(mi_array);

end