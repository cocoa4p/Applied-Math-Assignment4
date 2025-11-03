function dayTwo_Local()
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
    
    
    subplot(2,1,1);
    hold on;
    plot(t_range, V_list(:, 1), 'k', 'linewidth',2);
    plot(t_range, V_list(:, 2), 'b', 'linewidth', 2);

    plot(t_list, X_list(:,1),'r--', 'linewidth', 2);
    plot(t_list, X_list(:,2),'r--', 'linewidth', 2);
    xlabel('time');
    ylabel('position component');

    subplot(2,1,2);
    hold on;
    plot(t_range, V_list(:, 3), 'k', 'linewidth',2);
    plot(t_range, V_list(:, 4), 'b', 'linewidth', 2);
    
    plot(t_list, X_list(:,3),'r--', 'linewidth', 2);
    plot(t_list, X_list(:,4),'r--', 'linewidth', 2);
    xlabel('time');
    ylabel('velocity component');

    n_samples = 60;
    h_ref_list = logspace(-3, 1, n_samples);
    abs_diff_list = zeros(1,n_samples);

    tr_error_list = zeros(1,n_samples);
%     tr_error_list1 = zeros(1,n_samples);
%     tr_error_list2 = zeros(1,n_samples);

    for n = 1:length(h_ref_list)
        h_ref = h_ref_list(n);
        V_list = compute_planetary_motion(tspan(1)+h_ref,V0,orbit_params);
        
        % EXPLICIT RK STEP
        % Previously RK_step_embedded (rate_func_in,t,XA,h,BT_struct)

        [XB, num_evals] = explicit_RK_step(my_rate, tspan(1), V0, h_ref, expMethod);
        
        abs_diff_list(n) = norm(V_list - V0);
        tr_error_list(n) = norm(XB - V_list);

    
    end

    filter_params = struct();
    filter_params.min_yval = 1e-12;
    filter_params.max_yval = 1e-6;

    [p1,k1] = loglog_fit(h_ref_list, tr_error_list, filter_params);
   
    p1

    figure(2);

    loglog(h_ref_list, abs_diff_list, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 4);
    hold on
    loglog(h_ref_list, tr_error_list, 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 4);

    loglog(h_ref_list, k1*h_ref_list.^p1, 'm', 'LineWidth',1.5);

    title('Scaling of Local Truncation Error with Step Size')
    xlabel('Step size, h')
    ylabel('Error magnitude')
    legend('|X(t+h) - X(t)| true change','Local truncation error', 'Fit Line');


end