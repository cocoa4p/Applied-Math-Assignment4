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

function [t_list,X_list,h_avg, num_evals] = explicit_RK_fixed_step_integration(rate_func_in,tspan,X0,h_ref,BT_struct)
    N = ceil((tspan(2) - tspan(1))/h_ref); % number of time steps
    h_avg = ((tspan(2) - tspan(1))/N);

    t_list = linspace(tspan(1),tspan(2), N+1)';
    X_list = zeros(N+1, length(X0)); % every row corresponds to our state s, col nums will be heigh of X

    X_list(1,:) = X0';
    num_evals = 0;

    for m = 1:N
         % iterate through our N time steps
         % at each time step, we compute value at next time step
        [XB, num_evals_temp] = explicit_RK_step(rate_func_in, t_list(m),X_list(m,:)', h_avg, BT_struct);
        num_evals = num_evals + num_evals_temp;
        X_list(m+1, :) = XB';

    end

end