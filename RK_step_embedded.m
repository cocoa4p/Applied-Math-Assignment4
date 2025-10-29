%This function computes the value of X at the next time step
%for any arbitrary embedded RK method
%INPUTS:
%rate_func_in: the function used to compute dXdt. rate_func_in will
% have the form: dXdt = rate_func_in(t,X) (t is before X)
%t: the value of time at the current step
%XA: the value of X(t)
%h: the time increment for a single step i.e. delta_t = t_{n+1} - t_{n}
%BT_struct: a struct that contains the Butcher tableau
% BT_struct.A: matrix of a_{ij} values
% BT_struct.B: vector of b_i values
% BT_struct.C: vector of c_i values
%OUTPUTS:
%XB1: the approximate value for X(t+h) using the first row of the Tableau
%XB2: the approximate value for X(t+h) using the second row of the Tableau
%num_evals: A count of the number of times that you called
% rate_func_in when computing the next step
function [XB1, XB2, num_evals] = RK_step_embedded(rate_func_in,t,XA,h,BT_struct)
<<<<<<< Updated upstream

=======
        
    % The tableau
>>>>>>> Stashed changes
    A = BT_struct.A;
    B = BT_struct.B;
    C = BT_struct.C;


    s = length(C); % number of stages
    K = zeros(length(XA), s); % store all k_i values

    for i = 1:s
 
        sum_val1 = K * A(i,:)'; % Computes the weighted sum
<<<<<<< Updated upstream
        K(:, i) = rate_func_in(t + C(i)*h, XA + h*sum_val1); % rate function
    end

    % Combine all k_i to get X_n+1
    XB1 = XA + h * (K * B(1,:)'); % uses B vector as they are the slope values
    XB2 = XA + h * (K * B(2,:)'); % uses B vector as they are the slope values
    num_evals = s;
end

=======
        

        K(:, i) = rate_func_in(t + C(i)*h, XA + h*sum_val1); % rate function
    end


    XB = XA + h * (K * B(1,:)'); % uses B vector as they are the slope values
    num_evals = s;

end
>>>>>>> Stashed changes
