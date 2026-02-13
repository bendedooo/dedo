function omega = omega_D_A_D(t)
%omega_D/A|D as given for Problem 7. Will serve as the input function for
%the ode45 solver.


omega = [cos(2*t);
    cos(2*t);
    0.025*t];

end