function [axp, ayp] = perturb_hook(x, y, xd, yd, mu, eps)
% perturb_hook  Example small perturbation hook.
%   [axp, ayp] = perturb_hook(x, y, xd, yd, mu, eps)
% Default example: small central 1/r^2 perturbation centered on the larger primary (-mu).
% Replace this function with your desired perturbation (radiation pressure, oblateness, thrust, etc).

    r1vec = [x + mu; y];
    r1 = norm(r1vec);
    if r1 == 0
        axp = 0; ayp = 0; return;
    end

    accel_mag = eps / r1^2;
    axp = -accel_mag * (r1vec(1)/r1);
    ayp = -accel_mag * (r1vec(2)/r1);
end
