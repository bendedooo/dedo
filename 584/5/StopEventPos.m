%% Event for finishing the simulation once the position is near a given value.
function [value, isterminal, direction] = StopEventPos(t, x)
    value      = x(6)-0.5;% round(norm([x(1)-x(3) x(2)-x(4)]));
    isterminal = 1;   % Stop the integration
    direction  = 0;
end