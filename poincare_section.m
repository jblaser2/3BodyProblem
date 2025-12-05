function [xcross, xdotcross] = poincare_section(Y)
% poincare_section  Find y=0 crossings with ydot > 0 (linear interpolation)
%   [xcross, xdotcross] = poincare_section(Y)
%   Y is Nx4 array: [x, y, xdot, ydot]
%
%   returns arrays of x and xdot at the crossing points.

    x = Y(:,1); y = Y(:,2); xd = Y(:,3); yd = Y(:,4);
    xcross = [];
    xdotcross = [];

    for k=1:length(y)-1
        if (y(k) <= 0 && y(k+1) > 0) || (y(k) < 0 && y(k+1) >= 0) || (y(k) >= 0 && y(k+1) < 0)
            % linear interpolation to find y=0 between k and k+1
            denom = (y(k+1) - y(k));
            if denom == 0
                continue;
            end
            tfrac = -y(k) / denom;
            x_at_cross = x(k) + tfrac*(x(k+1)-x(k));
            xd_at_cross = xd(k) + tfrac*(xd(k+1)-xd(k));
            yd_at_cross = yd(k) + tfrac*(yd(k+1)-yd(k));

            if yd_at_cross > 0
                xcross(end+1,1) = x_at_cross; %#ok<AGROW>
                xdotcross(end+1,1) = xd_at_cross; %#ok<AGROW>
            end
        end
    end
end
