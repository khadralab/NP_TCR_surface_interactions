function y = plot_multivariate_MI(ax, activation, fname)
    surfaces = [1,3,5,10,20];
    valence_vals = [1,2,3,4,5,10];
    rho_vals = linspace(-4,3.5,16);
    rho_vals = 10.^ rho_vals;
    koff_vals = [0.0001, 0.001, 0.01, 0.1];
    np_radius = 20;
    threshold_array = [0,15,30,100];
    
    
    PI = 1/(length(valence_vals) * length(rho_vals) * length(koff_vals));
    input_entropy = log2(1/PI);
    mi = zeros(length(threshold_array),length(surfaces));
    
    if exist(fname, "file")
        load(fname);
    else
        for th = 1:length(threshold_array)
            threshold = threshold_array(th);
            for t = 1:length(surfaces)
                TPC = surfaces(t);
                joint_H = 0;
                input_vector = [];
                output_vector = [];
            
            
                for i = 1:length(valence_vals)
                    np_valence = valence_vals(i);
            
                    for j = 1:length(rho_vals)
                        np_rho = rho_vals(j);
            
                        for koff = koff_vals
                            [bound_tcr, bound_np] = load_BoundTCR_dist(koff, np_rho, np_radius, np_valence, TPC);
                            
                            if threshold == 0
                                digital_output = bound_tcr;
                            else
                                digital_output = digital_activation(bound_tcr, threshold);
                            end
            
                            nbins = max(digital_output) - min(digital_output) + 1;
                            counts = histcounts(digital_output, nbins);
                            
                            % Conditional probability of output O given input I
                            POgI = counts ./ sum(counts);
                            POgI = POgI(POgI~=0);
            
                            % Joint probability of output O and input I
                            POaI = POgI * PI;
                            
                            joint_H = joint_H - sum(POaI.*log2(POaI));
            
                            output_vector = [output_vector; digital_output];
            
                        end
                    end
                end
            
                output_entropy = shannon_entropy(output_vector);
                mi(th,t) = output_entropy + input_entropy - joint_H;
            end
        end
        save(fname,"mi");
    end

    disp(string(threshold_array(activation)));
    plot(ax, surfaces, mi(activation,:),':x','LineWidth',2); %,'DisplayName',);
    grid(ax,'on');
    xlabel(ax,'TPC');
    ylabel(ax,'MI');
    l = legend(ax, string(threshold_array(activation)));
    title(ax, "Threshold");
end