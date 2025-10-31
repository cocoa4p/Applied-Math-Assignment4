function [XB, num_evals] = explicit_RK_step(rate_func_in, t, XA, h, BT_struct)
 
    A = BT_struct.A;
    B = BT_struct.B(:); 
    C = BT_struct.C;

    s = length(C); % number of stages
    K = zeros(length(XA), s); % store all k_i values

    for i = 1:s
        % Weighted sum of previous K terms (explicit RK)
        sum_val1 = K(:,1:i-1) * A(i,1:i-1)';  
        % Evaluate rate function
        K(:,i) = rate_func_in(t + C(i)*h, XA + h*sum_val1);
    end

    % Combine all to get next X value
    XB = XA + h * (K * B);
    num_evals = s;
end
