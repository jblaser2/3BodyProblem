function [value,isterminal,direction] = yCrossEvent(~, y)
% yCrossEvent   ODE event function used to detect y = 0 crossings
% returns value = y. Not terminal; direction = 0 (both directions).
    value = y(2);        % y
    isterminal = 0;      % do not stop integration
    direction = 0;       % detect crossings in both directions
end
