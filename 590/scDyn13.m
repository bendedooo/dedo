function [out, magField, T_SC_c, magDipMom] = scDyn13(t,X,TimeStamp0,magFieldFlag,controlFlag,constRot)
% farkli orbitlerde eci to p yapilmasi icin orbit parameterlar hesaplanmali

% constants
Re  = 6378.2;
mu  = 398600.4405;
f = 0.00335; % flattening constant
J = [27 0 0; 0 17 0; 0 0 25];%./30; %kgm2

% get the states
r  = X(1:3);          % position vector (km)
rn = norm(r);
v  = X(4:6);          % velocity vector (km/sec)
omega = X(7:9);       % angular velocity vector (rad/sec)
q = X(10:13);         % qaternion:[q_v q_4]

% for magnetic field

if magFieldFlag == 1

    magField = getMagFieldNED(t,r,TimeStamp0); %[XYZ_NED' XYZ_ECEF XYZ_ECI XYZ_P];

    if controlFlag == 0

        T_mag = [0.01; 0; 0];
        T_SC_c = T_mag;

    elseif controlFlag ~= 0
        % Lovera's controller
        b_I = magField(:,3)*10^-9;
        O_I2B = qToDCM(q);
        b_B = O_I2B * b_I;
        b = b_B./norm(b_B);
        Gamma = eye(3)-b*b';

        % cotroller parameters
        Eps = 0.001;
        k_p = 10;
        k_v = 10;
        beta = 0.15;

        if controlFlag == 1 % saturated
            u = - (Eps^2*k_p*q(1:3) + Eps * beta * (min(1, max(-1, k_v*J*omega/beta))) );
        elseif controlFlag == 2 % unsaturated
            u = - (Eps^2*k_p*q(1:3) + Eps*k_v*J*omega );
        end
        T_SC_c =  Gamma * u;
        magDipMom = ssym(b_I)*u./(norm(b_I))^2;

    end
else
    magField = 0;
    T_mag = [0.01; 0; 0];
    T_SC_c = T_mag;
    magDipMom = [0;0;0];
end


% forces and torques on the system
Fpot   =  -mu/rn^3*r;
% Fpert  = [0 0 0]; %(1e-05)*e_v;
T_ext  = 0;

% state dots
if constRot ==1
    omega_dot = zeros(1,3);
else
    omega_dot = (J)\((-ssym(omega)*(J*omega)+T_SC_c+T_ext));
end
q_dot = 0.5.*ssym4(omega)*q;

Xdot(1:3)    = v;
Xdot(4:6)    = Fpot;
Xdot(7:9)    = omega_dot;
Xdot(10:13)  = q_dot;

out          = Xdot(:); % convert to vector column format

end