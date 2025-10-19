%This function computes the value of X at the next time step
%for any arbitrary RK method
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
%XB: the approximate value for X(t+h) (the next step)
% formula depends on the integration method used
%num_evals: A count of the number of times that you called
% rate_func_in when computing the next step

%   NOTES:
    % sum_val1 = K*(A(i,:)');
    % sum_val2 = K*B;

    % s = number of stages (how many k_i values you calculate
    % b_i = weights telling how much each slope k_i contributes to the final update

function [XB, num_evals] = explicit_RK_step(rate_func_in,t,XA,h,BT_struct)
    
    % The tableau
    A = BT_struct.A;
    B = BT_struct.B;
    C = BT_struct.C;


    s = length(C); % number of stages
    K = zeros(length(XA), s); % store all k_i values

    for i = 1:s
 
        sum_val1 = K * A(i,:)'; % Computes the weighted sum
        
        K(:, i) = rate_func_in(t + C(i)*h, XA + h*sum_val1); % rate function
    end

%     X_mid = XA + h/2*rate_func_in(t, XA); % midpoint estimate  
%     XB  = XA +  h*rate_func_in(t  + h/2, X_mid); % Uses midpoints slope 
  
    % Combine all k_i to get X_n+1
    XB = XA + h * (K * B'); % uses B vector as they are the slope values
    num_evals = s;

end
