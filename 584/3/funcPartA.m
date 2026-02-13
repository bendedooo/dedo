function x_out = funcPartA(t,x_in,u_func)
%3-2-1 Euler angles from equations given in class
phi = x_in(1);
theta = x_in(2);
psi = x_in(3);
u_in = u_func(t);

x_out = [1 sin(phi)*tan(theta) cos(phi)*tan(theta);
    0 cos(phi) -sin(phi)
    0 sin(phi)/cos(theta) cos(phi)/cos(theta)]*u_in;



end