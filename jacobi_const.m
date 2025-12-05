function C = jacobi_const(x,y,xd,yd,mu)
% jacobi_const  Compute Jacobi constant (vectorized).
%   C = jacobi_const(x,y,xd,yd,mu)
%
%   x,y,xd,yd may be scalars or column vectors of equal length.

    r1 = sqrt((x + mu).^2 + y.^2);
    r2 = sqrt((x - 1 + mu).^2 + y.^2);
    Omega = 0.5*(x.^2 + y.^2) + (1-mu)./r1 + mu./r2;
    C = 2.*Omega - (xd.^2 + yd.^2);
end
