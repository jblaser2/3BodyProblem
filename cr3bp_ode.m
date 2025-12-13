function dydt = cr3bp_ode(~, y, mu, perturb_handle, eps)
% cr3bp_ode  Equations of motion for the planar CR3BP in rotating frame.
%   dydt = cr3bp_ode(t, y, mu, perturb_handle, eps)
%   y = [x; y; xdot; ydot]
%
%   perturb_handle is a function handle: [axp, ayp] = perturb_handle(x,y,xd,yd,mu,eps)
%   if no perturbation is desired, set eps = 0 or pass [] for perturb_handle.

    x  = y(1); yy = y(2); u = y(3); v = y(4);
    % distances to primaries:
    r1 = sqrt((x + mu)^2 + yy^2);
    r2 = sqrt((x - 1 + mu)^2 + yy^2);

    % effective potential derivatives (Omega_x, Omega_y)
    Omega_x = x - (1-mu)*(x+mu)/r1^3 - mu*(x-1+mu)/r2^3;
    Omega_y = yy - (1-mu)*yy/r1^3 - mu*yy/r2^3;

    % Equations of motion in rotating frame:
    xdd =  2*v + Omega_x; 
    ydd = -2*u + Omega_y;

    % apply optional small perturbation
    if nargin >= 5 && eps ~= 0 && ~isempty(perturb_handle)
        [axp, ayp] = perturb_handle(x, yy, u, v, mu, eps);
        xdd = xdd + axp;
        ydd = ydd + ayp;
    end

    dydt = [u; v; xdd; ydd];
end
