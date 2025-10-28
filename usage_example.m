
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
    DormandPrince.B = [5179/57600, 0, 7571/16695, 393/640, -92097/339200, 187/2100, 1/40];
    DormandPrince.A = [0,0,0,0,0,0,0;
                        1/5, 0, 0, 0,0,0,0;...
                        3/40, 9/40, 0, 0, 0, 0,0;...
                        44/45, -56/15, 32/9, 0, 0, 0,0;...
                        19372/6561, -25360/2187, 64448/6561, -212/729, 0, 0,0;...
                        9017/3168, -355/33, 46732/5247, 49/176, -5103/18656, 0,0;...
                        35/384, 0, 500/1113, 125/192, -2187/6784, 11/84,0];

    h_ref = 0.01;
    [t_list, X_list, h_avg, num_evals] = explicit_RK_fixed_step_integration(my_rate, tspan, V0, h_ref, DormandPrince);
    
    
    figure(1);
    subplot(2,1,1);
    hold on;
    plot(t_range, V_list(:, 1), 'k', 'linewidth',2);
    plot(t_range, V_list(:, 2), 'b', 'linewidth', 2);

    plot(t_list, X_list(:,1),'r--', 'linewidth', 2);
    plot(t_list, X_list(:,2),'r--', 'linewidth', 2);


    xlabel('time');
    ylabel('position component');
    
    figure(2);
    subplot(2,1,1);
    hold on;
    plot(t_range, V_list(:, 3), 'k', 'linewidth',2);
    plot(t_range, V_list(:, 4), 'b', 'linewidth', 2);
    
    plot(t_list, X_list(:,3),'r--', 'linewidth', 2);
    plot(t_list, X_list(:,4),'r--', 'linewidth', 2);

    xlabel('time');
    ylabel('velocity component');
 
  
end