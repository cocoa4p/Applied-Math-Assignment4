function conservation_experiment()

    orbit_params = struct();
    orbit_params.m_sun = 1;
    orbit_params.m_planet = 1;
    orbit_params.G = 40;
    
    x0 = 8;
    y0 = 0;
    dxdt0 = 0;
    dydt0 = 1.5;

    X0 = [x0; y0; dxdt0;dydt0];

    E_ref = calc_mech_energy(orbit_params, X0);

    t0 = 0;
    tf = 10;
     
    rate = @(t, V) gravity_rate_func(t, V, orbit_params);

   % X0 = sol(t0);  
   % X_true = sol(tf);   

    forEu = @forward_euler_fixed_step_integration;

    h_step = logspace(-5, -0.1, 50); 

    errs_FE  = zeros(size(h_step));
    nfe_FE   = zeros(size(h_step));

    % Step Size Loop
    for r = 1:length(h_step)
        h = h_step(r);

        % Forward Euler
        [t_FE, X_FE, h_avg_FE, num_evals_FE] = forEu(rate, [t0 tf], X0, h);
        %X_FE_final = X_FE(end,:)';

        %errs_FE(r) = norm(X_FE_final - X_true);
        %nfe_FE(r) = num_evals_FE;

        % Energy tracking
        E_values = zeros(size(t_FE));
        for i = 1:length(t_FE)
            E_values(i) = calc_mech_energy(orbit_params, X_FE(i,:)');
        end
        E0 = E_values(1);
        rel_error = abs(E_values - E0) / abs(E0);

        % Record max relative energy deviation
        errs_FE(r) = max(rel_error);
        nfe_FE(r) = num_evals_FE;

    end

    % Plot energy conservation vs step size
    figure(1); clf
    loglog(h_step, errs_FE,  'ro', 'MarkerFaceColor','r', 'MarkerSize',2); hold on

    %filter_params = struct();
    %filter_params.max_xval = 1;

   % [p_FE, k_FE]   = loglog_fit(h_step, errs_FE, filter_params);
   % yline(E_ref)

    %loglog(h_step, k_FE*h_step.^p_FE, 'r-', 'LineWidth', 1);

    xlabel('step size h');
    ylabel('Max relative energy error');
    title('Energy Conservation vs Step Size (Forward Euler)');

    % Optionally fit slope on log-log plot
    filter_params = struct();
    filter_params.max_xval = 1;
    [p_FE, k_FE] = loglog_fit(h_step, errs_FE, filter_params);

    loglog(h_step, k_FE*h_step.^p_FE, 'r--', 'LineWidth', 1);
    legend(sprintf('Forward Euler (p = %.2f)', p_FE), 'Location','southwest');
end


function E_ref = calc_mech_energy(orbit_params, V)

    m_sun = orbit_params.m_sun;
    m_g = orbit_params.G;
    m_planet = orbit_params.m_planet;

    x = V(1);
    y = V(2);
    vx = V(3);
    vy = V(4);

    r = sqrt(x^2 + y^2);
    
    E_ref = 0.5 * m_sun  * (vx^2 + vy^2)  - m_g * m_sun * m_planet / r;

end