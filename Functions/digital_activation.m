function  y = digital_activation(bound_tcrs, threshold)
    y = max(0, bound_tcrs - threshold);
    %{
    if bound_tcrs > threshold
        y=1;
    else
        y=0;
    end
    %}
end