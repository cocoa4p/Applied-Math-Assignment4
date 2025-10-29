
%example for how to use compute_planetary_motion(...)
function usage_example()
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

    DormandPrince = struct();
    DormandPrince.C = [0, 1/5, 3/10, 4/5, 8/9, 1, 1];
    
    % DormandPrince.B = [5179/57600, 0, 7571/16695, 393/640, -92097/339200, 187/2100, 1/40];

    DormandPrince.B = [35/384, 0, 500/1113, 125/192, -2187/6784, 11/84, 0;...
5179/57600, 0, 7571/16695, 393/640, -92097/339200, 187/2100, 1/40];


    DormandPrince.A = [0,0,0,0,0,0,0;
                        1/5, 0, 0, 0,0,0,0;...
                        3/40, 9/40, 0, 0, 0, 0,0;...
                        44/45, -56/15, 32/9, 0, 0, 0,0;...
                        19372/6561, -25360/2187, 64448/6561, -212/729, 0, 0,0;...
                        9017/3168, -355/33, 46732/5247, 49/176, -5103/18656, 0,0;...
                        35/384, 0, 500/1113, 125/192, -2187/6784, 11/84,0];

    h_ref = 0.01;
    [t_list, X_list, h_avg, num_evals] = explicit_RK_fixed_step_integration(my_rate, tspan, V0, h_ref, DormandPrince);
    
  
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
 
    % LOCAL truncation error for embedded -------------------------------

    n_samples = 60;
    h_ref_list = logspace(-3, 1, n_samples);
    abs_diff_list = zeros(1,n_samples);
    tr_error_list1 = zeros(1,n_samples);
    tr_error_list2 = zeros(1,n_samples);

    for n = 1:length(h_ref_list)
        h_ref = h_ref_list(n);
        V_list = compute_planetary_motion(tspan(1)+h_ref,V0,orbit_params);
        
        [XB1, XB2, ~] = RK_step_embedded(my_rate, tspan(1), V0, h_ref, DormandPrince);
        
        abs_diff_list(n) = norm(V_list - V0);
        tr_error_list1(n) = norm(XB1 - V_list);
        tr_error_list2(n) = norm(XB2 - V_list);
    
    
    end

    filter_params = struct();
    filter_params.min_yval = 1e-12;
    filter_params.max_yval = 1e-6;

    [p1,k1] = loglog_fit(h_ref_list, tr_error_list1, filter_params);
    [p2,k2] = loglog_fit(h_ref_list, tr_error_list2, filter_params);

    p1
    p2

    figure(2);
    loglog(h_ref_list, abs_diff_list, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 2);
    hold on
    loglog(h_ref_list, tr_error_list1, 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 2);
    loglog(h_ref_list, tr_error_list2, 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 2);


    loglog(h_ref_list, k1*num_evals.^p1, 'k--', 'LineWidth',1.5);
    loglog(h_ref_list, k2*num_evals.^p2, 'b--', 'LineWidth',1.5);

%     GLOBAL truncation error -----------------------------------------

%     n_samples = 30;
%     h_ref_list = logspace(-3.3, 1, n_samples);
%     
%     num_evals_list = zeros(1,n_samples);
%     h_avg_list = zeros(1,n_samples);
%     tr_error_list = zeros(1, n_samples); % value of the global truncation error
% 
% 
%     for n = 1:length(h_ref_list)
%         
%         h_ref = h_ref_list(n);
% 
%         [t_list, X_list, h_avg, num_evals] = explicit_RK_fixed_step_integration(my_rate, tspan, V0, h_ref, DormandPrince);
%         tr_error = norm(X_list(end,:) - V_list(end,:));
%         tr_error_list(n) = tr_error;
%         h_avg_list(n) = h_avg;
%         num_evals_list(n) = num_evals;
%     end

%     filter_params = struct();
%     filter_params.min_yval = 1e-10;
%     filter_params.max_yval = 1;
%     
% 
%     [p1,k1] = loglog_fit(h_avg_list, tr_error_list, filter_params);
%     [p2,k2] = loglog_fit(num_evals_list, tr_error_list, filter_params);
%     
%     p1 = abs(p1);
%     p2 = abs(p2);
%     
%     figure(2);
%     loglog(h_avg_list, tr_error_list, 'ro', 'MarkerFaceColor','r');
%     hold on;
%     loglog(h_avg_list, k1*h_avg_list.^p1, 'r--', 'LineWidth',1.5);
%     
%     title('Local Truncation Error Plot as a Function of the Timestep')
%     xlabel('H average');
%     ylabel('Errors');
% 
%     figure(3);
%     loglog(num_evals_list, tr_error_list, 'bo', 'MarkerFaceColor','b');
%     hold on;
%     loglog(num_evals_list, k2*num_evals_list.^p2, 'b--', 'LineWidth',1.5);
%     title('Local Truncation Error Plot as a Function of the Number of Function Evaluations');
%     xlabel('Number of Evals');
%     ylabel('Errors');
%   
end