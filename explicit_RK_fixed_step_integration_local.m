%Runs numerical integration arbitrary RK method
%INPUTS:
%rate_func_in: the function used to compute dXdt. rate_func_in will
% have the form: dXdt = rate_func_in(t,X) (t is before X)
%tspan: a two element vector [t_start,t_end] that denotes the integration endpoints
%X0: the vector describing the initial conditions, X(t_start)
%h_ref: the desired value of the average step size (not the actual value)
%BT_struct: a struct that contains the Butcher tableau
% BT_struct.A: matrix of a_{ij} values
% BT_struct.B: vector of b_i values
% BT_struct.C: vector of c_i values
%OUTPUTS:
%t_list: the vector of times, [t_start;t_1;t_2;...;.t_end] that X is approximated at
%X_list: the vector of X, [X0';X1';X2';...;(X_end)'] at each time step
%h_avg: the average step size
%num_evals: total number of calls made to rate_func_in during the integration

function [t_list,X_list,h_avg, num_evals] = explicit_RK_fixed_step_integration_local(rate_func_in,tspan,X0,h_ref,BT_struct)
    % Initialize
    t = tspan(1);
    X = X0(:);
    t_list = t;
    X_list = X.';
    num_evals = 0;
    
    h = h_ref; % initial step size
    tol = 1e-6; % desired local error tolerance
    
    while t < tspan(2)
    if t + h > tspan(2)
    h = tspan(2) - t; % adjust last step to hit t_end exactly
    end

    % Take a single RK step
    [XB, num_evals_temp, err_est] = explicit_RK_step(rate_func_in, t, X, h, BT_struct); 
    num_evals = num_evals + num_evals_temp;
    
    % Compute new step size using simple error control
    if err_est <= tol
        % Accept step
        t = t + h;
        X = XB;
        t_list(end+1,1) = t;
        X_list(end+1,:) = X.';
    end
    
    % Update step size for next iteration (safety factor 0.9)
    h = 0.9 * h * (tol / max(err_est,1e-16))^(1/5); 

end


end