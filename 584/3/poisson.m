function x_dot = poisson(t,x_in,u_func)
x_mat = reshape(x_in,3,3);
u  = u_func(t);
x_dot_mat = -[0 -u(3) u(2);u(3) 0 -u(1);-u(2) u(1) 0]*x_mat;
x_dot = reshape(x_dot_mat,9,1);
end