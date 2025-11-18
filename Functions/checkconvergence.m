function [convergence] = checkconvergence(x, window, tau)
%x = array representing some time-series of state variable
%window = index interval over which mean and variance are estimated
%tau = time step between elements of array x

    window = floor(window / tau);
    
    if length(x) < 2*window+1
        convergence = false;
        return
    end
    
    x1 = x(end-2*window:end-window);
    x2 = x(end-window:end);
    delta_mean = abs(mean(x2) - mean(x1)) / mean(x2);
    
    delta_std = abs(std(x2) - std(x1)) / std(x2);
    
    if delta_mean < 0.05 & delta_std < 1
        convergence = true;
        return
    else
        convergence = false;
        return
    end
    
end