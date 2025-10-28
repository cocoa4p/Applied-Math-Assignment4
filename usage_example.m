
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
    
    V0 = [x0;y0;dxdt0;dydt0];
    t_range = linspace(0,30,100);
    V_list = compute_planetary_motion(t_range,V0,orbit_params);
    
    % WIP SEE VIDEO
%     hold on;
%     subplot(2,1,1);
%     plot(t_range, V_list(:, 1))
%     plot(t_range, V_list(:, 1))
 
    % WIP ABOVE ________


    axis equal; axis square;
    axis([-20,20,-20,20])
    hold on
    plot(0,0,'ro','markerfacecolor','r','markersize',5);
    plot(V_list(:,1),V_list(:,2),'k');
end