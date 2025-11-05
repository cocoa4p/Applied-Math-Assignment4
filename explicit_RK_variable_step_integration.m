%Runs numerical integration arbitrary RK method using variable time steps
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
%p: how error scales with step size (error = k*hˆp)
%error_desired: the desired local truncation error at each step
%OUTPUTS:
%t_list: the vector of times, [t_start;t_1;t_2;...;.t_end] that X is approximated at
%X_list: the vector of X, [X0';X1';X2';...;(X_end)'] at each time step
%h_avg: the average step size
%num_evals: total number of calls made to rate_func_in during the integration

function [t_list,X_list,h_avg, num_evals,fail_rate] = explicit_RK_variable_step_integration ...
(rate_func_in,tspan,X0,h_ref,BT_struct,p,error_desired)
    t = tspan(1);
    X = X0;
    t_end = tspan(2);

    t_list = t;
    X_list = X.';
    num_evals = 0;
    h = h_ref;
    n_attempts = 0;
    n_fail = 0;

    while t < t_end
        if t + h > t_end
            h = t_end - t;
        end
        [XB, n_eval, h_next, redo] = explicit_RK_variable_step( ...
            rate_func_in, t, X, h, BT_struct, p, error_desired);
        n_attempts = n_attempts + 1;
        num_evals = num_evals + n_eval;

        if ~redo
            t = t + h;
            X = XB;
            t_list(end+1,1) = t;
            X_list(end+1,:) = X.';
        else
            n_fail = n_fail + 1;
        end
        h = h_next;
    end
    h_avg = (t_list(end)-t_list(1))/(length(t_list)-1);
    fail_rate = n_fail / n_attempts;
    fprintf('Step failure rate = %.3f\n', fail_rate);
end
