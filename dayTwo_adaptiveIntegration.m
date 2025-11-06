function dayTwo_adaptiveIntegration()

    % ------------------ SETUP ------------------
    orbit_params = struct();
    orbit_params.m_sun = 1;
    orbit_params.m_planet = 1;
    orbit_params.G = 50;

    % Initial conditions (elliptical orbit)

    x0 = 10;
    y0 = 0;
    dxdt0 = 0;
    dydt0 = 3;   % smaller tangential velocity


%     x0 = 8;  y0 = 0;
%     dxdt0 = 0;  dydt0 = 1.5;
    V0 = [x0; y0; dxdt0; dydt0];

    tspan = [0, 400];
    t_range = linspace(tspan(1), tspan(2), 500);
    V_true = compute_planetary_motion(t_range, V0, orbit_params);

    % Rate function handle
    my_rate = @(t, V) gravity_rate_func(t, V, orbit_params);

    % Dormand–Prince Butcher tableau (embedded RK 5(4))
    % Embedded version (for adaptive step)
    DormandPrince_embedded = struct();
    DormandPrince_embedded.C = [0, 1/5, 3/10, 4/5, 8/9, 1, 1];
    DormandPrince_embedded.B = [ ...
        35/384, 0, 500/1113, 125/192, -2187/6784, 11/84, 0; ...  % higher-order (5th)
        5179/57600, 0, 7571/16695, 393/640, -92097/339200, ...
        187/2100, 1/40];                                          % lower-order (4th)
    DormandPrince_embedded.A = [ ...
        0,0,0,0,0,0,0; ...
        1/5, 0,0,0,0,0,0; ...
        3/40, 9/40, 0,0,0,0,0; ...
        44/45, -56/15, 32/9, 0,0,0,0; ...
        19372/6561, -25360/2187, 64448/6561, -212/729, 0,0,0; ...
        9017/3168, -355/33, 46732/5247, 49/176, -5103/18656, 0,0; ...
        35/384, 0, 500/1113, 125/192, -2187/6784, 11/84, 0];

    % Single-row version (for fixed step)
    DormandPrince_fixed = struct();
    DormandPrince_fixed.C = DormandPrince_embedded.C;
    DormandPrince_fixed.A = DormandPrince_embedded.A;
    DormandPrince_fixed.B = DormandPrince_embedded.B(2,:);  % choose one row
    
    

    p = 5; 


    error_desired_list = logspace(-3, -14, 50);  % adaptive target errors
    h_fixed_list = logspace(-2.5, 1, 50);      % fixed step reference sizes

    global_error_fixed = zeros(size(h_fixed_list));
    h_avg_fixed = zeros(size(h_fixed_list));
    num_eval_fixed = zeros(size(h_fixed_list));

    global_error_adapt = zeros(size(error_desired_list));
    h_avg_adapt = zeros(size(error_desired_list));
    num_eval_adapt = zeros(size(error_desired_list));
    fail_rate_list = zeros(size(error_desired_list));


    fprintf('Running fixed-step integrations...\n');
    for i = 1:length(h_fixed_list)
        h = h_fixed_list(i);
        [t_list, X_list, h_avg, num_evals] = explicit_RK_fixed_step_integration_global( ...
            my_rate, tspan, V0, h, DormandPrince_fixed);
        global_error_fixed(i) = norm(X_list(end,:) - V_true(end,:));
        h_avg_fixed(i) = h_avg;
        num_eval_fixed(i) = num_evals;
    end


    fprintf('Running adaptive-step integrations...\n');
    for i = 1:length(error_desired_list)
        error_desired = error_desired_list(i);
        [t_list, X_list, h_avg, num_evals, fail_rate] = explicit_RK_variable_step_integration(my_rate, tspan, V0, 0.05, DormandPrince_embedded, p, error_desired);
        global_error_adapt(i) = norm(X_list(end,:) - V_true(end,:));
        h_avg_adapt(i) = h_avg;
        num_eval_adapt(i) = num_evals;
        fail_rate_list(i) = fail_rate;
    end

    figure(1); clf;
    loglog(h_avg_fixed, global_error_fixed, 'bo-', 'LineWidth', 1.5, 'DisplayName', 'Fixed step');
    hold on;
    loglog(h_avg_adapt, global_error_adapt, 'ro-', 'LineWidth', 1.5, 'DisplayName', 'Adaptive step');
    xlabel('Average step size, h');
    ylabel('Global truncation error');
    legend('Location', 'best');
    title('Global Error vs. Average Step Size');

    
    figure(2); clf; 
    loglog(num_eval_fixed, global_error_fixed, 'bo-', 'LineWidth', 1.5, 'DisplayName', 'Fixed step');
    hold on;
    loglog(num_eval_adapt, global_error_adapt, 'ro-', 'LineWidth', 1.5, 'DisplayName', 'Adaptive step');
    xlabel('Number of function evaluations');
    ylabel('Global truncation error');
    legend('Location', 'best');
    title('Global Error vs. Number of Evaluations');


    figure(3); clf;
    
    [t_list, X_list, ~, ~] = explicit_RK_variable_step_integration( ...
        my_rate, tspan, V0, 0.05, DormandPrince_embedded, p, 1e-5);
    plot(X_list(:,1), X_list(:,2), 'r-', 'LineWidth', 1.5);
    hold on;
    plot(V_true(:,1), V_true(:,2), 'k--', 'LineWidth', 1.2);
    % xlim([-80 100])
    % ylim([-80 100])
    axis equal;
    plot(0,0,'ro','markerfacecolor','r','markersize',5);
    xlabel('x'); ylabel('y');
    legend('Adaptive RK (Dormand–Prince)', 'True solution');
    title('Planetary Orbit Comparison');
    axis equal;


    figure(4); clf;
    h_list = diff(t_list);
    r_list = sqrt(sum(X_list(1:end-1,1:2).^2, 2));
    semilogx(r_list, h_list, 'ko', 'MarkerFaceColor', 'b', 'MarkerSize', 4);

    figure(5); clf;
    semilogx(h_avg_adapt, fail_rate_list, 'co', 'MarkerFaceColor', 'c', 'MarkerSize', 2)
    xlabel('average step size (log)'); ylabel('failure rate');
    %legend('Adaptive RK (Dormand–Prince)', 'True solution');
    title('Failure Rate Over Average Step Size');



    fprintf('Running adaptive-step integration...\n');
        error_desired = 1e-8;
        [t_list, X_list, h_avg, num_evals, fail_rate] = explicit_RK_variable_step_integration(my_rate, tspan, V0, 0.05, DormandPrince_embedded, p, error_desired);
        %global_error_adapt(i) = norm(X_list(end,:) - V_true(end,:));
        %h_avg_adapt(i) = h_avg;
        %num_eval_adapt(i) = num_evals;
        %fail_rate_list(i) = fail_rate;


    % Position vs time
    figure(6);
    plot(t_list, X_list(:,1), 'r-','marker', 'o', 'markerfacecolor', 'k', 'markeredgecolor', 'k', 'markersize', 2);
    xlabel('Time');
    ylabel('Position');
    title('Adaptive RK Integration: Position vs Time');
    
    % Velocity vs time
    figure(7);
    plot(t_list, X_list(:,2), 'r-', 'marker', 'o', ...
         'markerfacecolor', 'k', ...
         'markeredgecolor', 'k', ...
         'markersize', 2);
    xlabel('Time');
    ylabel('Velocity');
    title('Adaptive RK Integration: Velocity vs Time');


    % step sizes between adjacent times
    h_list = diff(t_list);
    
    % Compute distances from the Sun (assuming X_list(:,1)=x, X_list(:,2)=y)
    r_list = sqrt(X_list(1:end-1,1).^2 + X_list(1:end-1,2).^2);
    
    % Scatter plot of step size vs distance
    figure(8);
    scatter(r_list, h_list, 10, 'filled');
    xlabel('Distance from Sun, r');
    ylabel('Step size, h');
    title('Adaptive RK: Step size vs Distance');
    grid on;