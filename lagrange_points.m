function [Lx, Ly] = lagrange_points(mu)
% lagrange_points  Compute locations of the five classical Lagrange points
%   [Lx, Ly] = lagrange_points(mu)
%   returns arrays Lx(1..5), Ly(1..5) corresponding to L1..L5.

    % function for x-derivative when y=0
    fx = @(x) x - (1-mu)*(x+mu)./abs(x+mu).^3 - mu*(x-1+mu)./abs(x-1+mu).^3;

    % initial guesses chosen according to standard CR3BP geometry
    L1x = fzero(fx, 0.7);
    L2x = fzero(fx, 1.2);
    L3x = fzero(fx, -1.0);

    % triangular points (equilateral triangle positions)
    xL4 = 0.5 - mu;
    yL4 =  sqrt(3)/2;
    xL5 = 0.5 - mu;
    yL5 = -sqrt(3)/2;

    Lx = [L1x; L2x; L3x; xL4; xL5];
    Ly = [0; 0; 0; yL4; yL5];
end
