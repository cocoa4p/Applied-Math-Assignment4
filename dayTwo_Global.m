function dayTwo_Global()
    orbit_params = struct();
    orbit_params.m_sun = 1;
    orbit_params.m_planet = 1;
    orbit_params.G = 40;
  
    x0 = 8;
    y0 = 0;
    dxdt0 = 0;
    dydt0 = 1.5;
    
    V0 = [x0; y0; dxdt0;dydt0];
    tspan = [0,30];
    t_range = linspace(tspan(1),tspan(2),100);
    V_list = compute_planetary_motion(t_range,V0,orbit_params);
    
    my_rate = @(t_in, V_in) gravity_rate_func(t_in, V_in, orbit_params);

    % Three different explicit Runge-Kutta Methods to test

    % Ralston's Method
    Ralston = struct();
    Ralston.C = [0; 2/3];
    Ralston.B = [1/4, 3/4];
    Ralston.A = [0, 0; 2/3, 0];

    % Heun's Third-Order Method
    Huen3 = struct();
    Huen3.C = [0; 0; 0];
    Huen3.B = [0; 1/3; 2/3];
    Huen3.A = [0, 0, 0; 1/3, 0, 0; 0, 2/3, 0];

    % "Original" fourth-order Runge-Kutta
    ogRunge = struct();
    ogRunge.C = [0; 1/2; 1/2; 1];
    ogRunge.B = [1/6, 1/3, 1/3, 1/6];
    ogRunge.A = [0, 0, 0, 0; 1/2, 0, 0, 0; 0, 1/2, 0, 0; 0, 0, 1, 0];

    
    DormandPrince = struct();
    DormandPrince.C = [0, 1/5, 3/10, 4/5, 8/9, 1, 1];
    % DormandPrince.B = [35/384, 0, 500/1113, 125/192, -2187/6784, 11/84, 0; 5179/57600, 0, 7571/16695, 393/640, -92097/339200, 187/2100, 1/40];
    DormandPrince.B = [5179/57600, 0, 7571/16695, 393/640, -92097/339200, 187/2100, 1/40];
    DormandPrince.A = [0,0,0,0,0,0,0;
    1/5, 0, 0, 0,0,0,0;...
    3/40, 9/40, 0, 0, 0, 0,0;...
    44/45, -56/15, 32/9, 0, 0, 0,0;...
    19372/6561, -25360/2187, 64448/6561, -212/729, 0, 0,0;...
    9017/3168, -355/33, 46732/5247, 49/176, -5103/18656, 0,0;...
    35/384, 0, 500/1113, 125/192, -2187/6784, 11/84,0];

% Global trucnation for the conserved. take the estimated at the end to the
% initial value to find the amount it changed. then run integrator at
% different time step sizes


%     expMethod = Ralston;
%     expMethod = Huen3;
%     expMethod = ogRunge;
    expMethod = DormandPrince;

    h_ref = 0.05;
    [t_list, X_list, h_avg, num_evals] = explicit_RK_fixed_step_integration(my_rate, tspan, V0, h_ref, expMethod);
    

% GLOBAL ERROR -------------------------
    n_samples = 30;
    h_ref_list = logspace(-3.3, 1, n_samples);
    
    num_evals_list = zeros(1,n_samples);
    h_avg_list = zeros(1,n_samples);
    tr_error_list = zeros(1, n_samples); % value of the global truncation error


    for n = 1:length(h_ref_list)
        
        h_ref = h_ref_list(n);

        [t_list, X_list, h_avg, num_evals] = explicit_RK_fixed_step_integration(my_rate, tspan, V0, h_ref, expMethod);
        tr_error = norm(X_list(end,:) - V_list(end,:));
        tr_error_list(n) = tr_error;
        h_avg_list(n) = h_avg;
        num_evals_list(n) = num_evals;
    end

    filter_params = struct();
    filter_params.min_yval = 1e-10;
    filter_params.max_yval = 1;
    

    [p1,k1] = loglog_fit(h_avg_list, tr_error_list, filter_params);
    [p2,k2] = loglog_fit(num_evals_list, tr_error_list, filter_params);
    
    p1 = abs(p1);
    p2 = abs(p2);
    
    figure(2);
    loglog(h_avg_list, tr_error_list, 'ro', 'MarkerFaceColor','r');
    hold on;
    loglog(h_avg_list, k1*h_avg_list.^p1, 'r--', 'LineWidth',1.5);
    
    title('Global Truncation Error vs Step Size');
    xlabel('Average timestep (h)');
    ylabel('Error');
    legend('Computed error', sprintf('Fit: slope = %.2f', p1), 'Location', 'best');

    figure(3);
    loglog(num_evals_list, tr_error_list, 'bo', 'MarkerFaceColor','b');
    hold on;
    loglog(num_evals_list, k2*num_evals_list.^p2, 'b--', 'LineWidth',1.5);
    title('Global Truncation Error vs Number of Evaluations');
    xlabel('Number of function evaluations');
    ylabel('Error');
    legend('Computed error', sprintf('Fit: slope = %.2f', p2), 'Location', 'best');


end