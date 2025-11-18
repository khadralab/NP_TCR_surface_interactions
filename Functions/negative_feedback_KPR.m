function [Cn, S, phos_coef, b, r2] = negative_feedback_KPR(bound_tcrs, koff, kp, ST, Cs, feedbacktrue)
    % Parameters
    N = kp(3);
    b = kp(2);
    kp = kp(1);
    %ST = 0.5;
    %Cs = 100e0;  % Threshold for 50% activation of phosphotase (number of TCRs)

    if feedbacktrue
        phosfun = @(S)phosphotase_solver(S, kp, koff, b, ST, Cs, bound_tcrs);
    
        s0 = 10;
    
        options  = optimoptions('fsolve','Display','off','OptimalityTolerance',1e-8,'StepTolerance',1e-8,...
            'FunctionTolerance',1e-9,'Algorithm','trust-region-dogleg','MaxFunctionEvaluations',200);
        [S, fval, exitflag, output] = fsolve(phosfun, s0, options)
    else
        S=0;
    end

    r1 = r_neg(S, kp, koff, b);
    r2 = r_pos(S, kp, koff, b);
    phos_coef = (1 - r1 / r2) * r1 ^ N;

    Cn = bound_tcrs * (1 - r1 / r2) * r1 ^ N;
%{
    ST = linspace(0, 5, 1000);
    S_array = zeros(size(ST));

    for i = 1:length(ST)
        phosfun = @(S)phosphotase_solver(S, kp, koff, b, ST(i), Cs, bound_tcrs);

        s0 = 300;
        S_array(i) = fsolve(phosfun, s0);
    end
    r1 = r_neg(S_array, kp, koff, b);
    r2 = r_pos(S_array, kp, koff, b);

    phos_coef = (1 - r1 / r2) .* r1 .^ N;

    Cn = bound_tcrs * (1 - r1 / r2) .* r1 .^ N;
    
    
    figure()
    for i=1:5
        ax(i) = subplot(3,2,i);
    end
    
    x=ST;
    plot(ax(1), x, Cn);
    plot(ax(2), x, phos_coef);
    plot(ax(3), x, r1);
    plot(ax(4), x, r2);
    plot(ax(5), x, S_array);

    title(ax(1), 'Phos. Receptors')
    title(ax(2), 'Phos Coef')
    title(ax(3), 'R negative')
    title(ax(4), 'R positive')
    title(ax(5), 'S')
    %}
end

function y = phosphotase_solver(S, kp, koff, b, ST, Cs, bound_tcrs)
    S = S/100;
    r = r_neg(S, kp, koff, b);
    y = ST.*r.*(1-r) ./ (r.*(1-r) + Cs ./ bound_tcrs) - S;
end

function y = r_neg(S, kp, koff, b)
    y = (kp + b + S + koff - discriminant(S,kp,koff,b)) ./ (2*b + 2*S);
end

function y = r_pos(S, kp, koff, b)
    y = (kp + b + S + koff + discriminant(S,kp,koff,b)) ./ (2*b + 2*S);
end

function y = discriminant(S, kp, koff, b)
    y = sqrt((kp + b + S + koff).^2 - 4 * kp * (b + S));
end