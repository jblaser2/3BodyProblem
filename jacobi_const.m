function C = jacobi_const(x,y,xd,yd,mu)
% jacobi_const  Compute Jacobi constant (vectorized).
%   C = jacobi_const(x,y,xd,yd,mu)
%
%   x,y,xd,yd may be scalars or column vectors of equal length.


% Jacobi const is sort of like energy conservation. It is constant along
% certain trajectories. A higher C means less kinetic energy available. 

    r1 = sqrt((x + mu).^2 + y.^2); % The distance to the Earth
    r2 = sqrt((x - 1 + mu).^2 + y.^2); % Distance to the moon
    Omega = 0.5*(x.^2 + y.^2) + (1-mu)./r1 + mu./r2; % The effective potential
    C = 2.*Omega - (xd.^2 + yd.^2); % The Jacobi constant- denoting the zero-velocity curves where 2.*Omega = C
end
