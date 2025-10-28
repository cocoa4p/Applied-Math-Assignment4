%Rate function describing Newton's law of gravitation
%INPUTS:
%t: the current time
%V: the vector of the position and velocity of the planet
% V = [x_p; y_p; dxdt_p; dydt_p]

%orbit_params: a struct describing the system parameters
% orbit_params.m_sun is the mass of the sun
% orbit_params.m_planet is the mass of the planet
% orbit_params.G is the gravitational constant
%OUTPUTS:
%dVdt: a column vector describing the time derivative of V:
% dVdt = [dxdt_p; dydt_p; d2xdt2_p; d2ydt2_p]
function dVdt = gravity_rate_func(t,V,orbit_params)
    
    x = V(1); y = V(2); vx = V(3); vy = V(4);
    r = [x;y];

%   r_mag = sqrt(x^2 + y^2);
    r_mag = norm(r);
    m_sun = orbit_params.m_sun;
    G = orbit_params.G;

    a = (-(m_sun*G)/r_mag^3)*r;
    ax = a(1); ay = a(2);

    dVdt = [vx; vy; ax; ay];


end
