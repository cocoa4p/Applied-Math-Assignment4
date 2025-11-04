function dayTwo_adaptiveIntegrationREVISED()

    orbit_params = struct('m_sun',1,'m_planet',1,'G',30);
    
    % Initial conditions for eccentric elliptical orbit
    x0 = 10; y0 = 0; dxdt0 = 0; dydt0 = 1.2;
    V0 = [x0; y0; dxdt0; dydt0];
    
    tspan = [0,30];
    t_ref = linspace(tspan(1), tspan(2), 500);
    V_true = compute_planetary_motion(t_ref, V0, orbit_params);  % reference solution
    
    % Rate function
    my_rate = @(t,V) gravity_rate_func(t,V,orbit_params);
    
    % Dormand-Prince embedded RK (5th/4th order)
    DP_embedded = getDormandPrinceEmbedded();
    
    % Fixed-step uses lower-order row of B
    DP_fixed = DP_embedded;
    DP_fixed.B = DP_embedded.B(2,:);  % lower order
    
    % Order of method
    p = 5;
    
    %% ------------------ LTE AND STEP SIZES ------------------
    h_list = logspace(-2.5,-1,5);            % candidate fixed step sizes
    error_desired_list = logspace(-3,-8,6);  % adaptive targets
    
    % Storage
    LTE_XB1 = cell(length(h_list),1);
    LTE_XB2 = cell(length(h_list),1);
    LTE_diff = cell(length(h_list),1);
    h_actual_list = cell(length(h_list),1);
    
    % Loop over fixed steps to compute local truncation errors
    for i = 1:length(h_list)
        h = h_list(i);
        [t_list, X_list, XB1_list, XB2_list, h_actual] = explicit_RK_fixed_step_integration_LTE(my_rate, tspan, V0, h, DP_embedded);
        h_actual_list{i} = h_actual;
    
        % Local truncation errors per step
        LTE_XB1{i} = vecnorm(XB1_list - V_true(1:length(XB1_list),:),2,2);
        LTE_XB2{i} = vecnorm(XB2_list - V_true(1:length(XB2_list),:),2,2);
        LTE_diff{i} = vecnorm(XB1_list - XB2_list,2,2);
    end
    
    %% ------------------ Adaptive RK with step failure rate ------------------
    n_adapt = length(error_desired_list);
    global_error_adapt = zeros(n_adapt,1);
    h_avg_adapt = zeros(n_adapt,1);
    num_eval_adapt = zeros(n_adapt,1);
    fail_rate_adapt = zeros(n_adapt,1);
    
    for i = 1:n_adapt
        error_desired = error_desired_list(i);
        [t_list, X_list, h_avg, num_evals, num_fail] = explicit_RK_variable_step_integration(my_rate, tspan, V0, 0.05, DP_embedded, p, error_desired);
        global_error_adapt(i) = norm(X_list(end,:) - V_true(end,:));
        h_avg_adapt(i) = h_avg;
        num_eval_adapt(i) = num_evals;
        fail_rate_adapt(i) = num_fail/num_evals;
    end
    
    %% ------------------ FIXED STEP GLOBAL ERRORS ------------------
    global_error_fixed = zeros(length(h_list),1);
    h_avg_fixed = zeros(length(h_list),1);
    num_eval_fixed = zeros(length(h_list),1);
    
    for i = 1:length(h_list)
        h = h_list(i);
        [t_list, X_list, ~, ~, h_avg, num_evals] = explicit_RK_fixed_step_integration_global(my_rate, tspan, V0, h, DP_fixed);
        global_error_fixed(i) = norm(X_list(end,:) - V_true(end,:));
        h_avg_fixed(i) = h_avg;
        num_eval_fixed(i) = num_evals;
    end
    
    %% ------------------ PLOTS ------------------
    % --- Local truncation error vs h ---
    figure; clf; hold on;
    for i=1:length(h_list)
        loglog(h_list(i)*ones(size(LTE_XB1{i})), LTE_XB1{i}, 'bo');
        loglog(h_list(i)*ones(size(LTE_XB2{i})), LTE_XB2{i}, 'rs');
        loglog(h_list(i)*ones(size(LTE_diff{i})), LTE_diff{i}, 'k.');
    end
    xlabel('Step size h'); ylabel('Local truncation error');
    legend('XB1','XB2','|XB1-XB2|'); title('LTE vs Step Size');
    
    % --- LTE vs |XB1-XB2| ---
    figure; clf; hold on;
    for i=1:length(h_list)
        loglog(LTE_diff{i}, LTE_XB1{i}, 'bo');
        loglog(LTE_diff{i}, LTE_XB2{i}, 'rs');
    end
    xlabel('|XB1-XB2|'); ylabel('Local truncation error');
    legend('XB1','XB2'); title('LTE vs Embedded RK difference');
    
    % --- Global error vs average step size ---
    figure; clf; hold on;
    loglog(h_avg_fixed, global_error_fixed,'bo-');
    loglog(h_avg_adapt, global_error_adapt,'rs-');
    xlabel('Average step size'); ylabel('Global error');
    legend('Fixed step','Adaptive step'); title('Global error vs Avg step size');
    
    % --- Global error vs function evaluations ---
    figure; clf; hold on;
    loglog(num_eval_fixed, global_error_fixed,'bo-');
    loglog(num_eval_adapt, global_error_adapt,'rs-');
    xlabel('Number of function evaluations'); ylabel('Global error');
    legend('Fixed step','Adaptive step'); title('Global error vs # evals');
    
    % --- Step failure rate vs average step ---
    figure; clf;
    semilogx(h_avg_adapt, fail_rate_adapt,'ro-','LineWidth',1.5);
    xlabel('Average step size'); ylabel('Step failure rate'); title('Adaptive step failure rate');
    
    % --- Position and velocity vs time for one adaptive run ---
    [t_list, X_list, ~, ~, ~] = explicit_RK_variable_step_integration(my_rate, tspan, V0, 0.05, DP_embedded, p, 1e-5);
    figure; clf; subplot(2,1,1);
    plot(t_list,X_list(:,1),'ro-','markerfacecolor','k','markersize',2); hold on;
    xlabel('Time'); ylabel('x-position'); title('Position vs time');
    subplot(2,1,2);
    plot(t_list,X_list(:,3),'bo-','markerfacecolor','k','markersize',2); hold on;
    xlabel('Time'); ylabel('x-velocity'); title('Velocity vs time');
    
    % --- Step-size clustering vs distance from star ---
    r_list = sqrt(sum(X_list(:,1:2).^2,2));
    h_list_adapt = diff(t_list);
    figure; clf;
    semilogx(r_list(1:end-1), h_list_adapt,'ko','markerfacecolor','b','markersize',4);
    xlabel('Distance from star r'); ylabel('Adaptive step size'); title('Step size vs radius');

end