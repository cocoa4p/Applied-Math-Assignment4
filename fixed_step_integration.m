%Runs fixed step numerical integration
%INPUTS:
%rate_func_in: the function used to compute dXdt. rate_func_in will
% have the form: dXdt = rate_func_in(t,X) (t is before X)
%step_func: the function used to approximate X(t) at the next step
% [XB,num_evals] = step_func(rate_func_in,t,XA,h)
%tspan: a two element vector [t_start,t_end] that denotes the integration endpoints
%X0: the vector describing the initial conditions, X(t_start)
%h_ref: the desired value of the average step size (not the actual value)
%OUTPUTS:
%t_list: the vector of times, [t_start;t_1;t_2;...;.t_end] that X is approximated at
%X_list: the vector of X, [X0';X1';X2';...;(X_end)'] at each time step
%h_avg: the average step size
%num_evals: total number of calls made to rate_func_in during the integration

function [t_list,X_list,h_avg, num_evals] = fixed_step_integration(rate_func_in,step_func,tspan,X0,h_ref)
    
    n_steps = ceil((tspan(2) - tspan(1))/h_ref);
    h_avg = (tspan(2) - tspan(1))/n_steps;
    
    t_list = linspace(tspan(1),tspan(2),n_steps+1)';
    X_list =  zeros(length(t_list),length(X0));
    
    Xval =  X0;
    X_list(1,:) = Xval;
    
    num_evals = 0;
    for n =  1:n_steps
        
        [Xval, num_evals_temp] = step_func(rate_func_in, t_list(n), Xval, h_avg);
        num_evals = num_evals + num_evals_temp;
     
        X_list(n+1,:) = Xval';
    end


end
