Re      =  6378;           % Earth Radius (km)
mu      =  398600.4;       % Gravitational constant for Earth (km^3/s^2)
r       =  Re + 350;       % Radius of the circular orbit
n       =  sqrt(mu/r^3);   % orbital rate / mean motion

% Continuous-time CW equation matrices
A = [0,     0,      1,      0;
     0,     0,      0,      1;
     3*n^2, 0,      0,      2*n;
     0,     0,      -2*n,   0]
 
 dT = 30; % period between two subsequent delta-v's
 
 % Discrete-time dynamics
 Ad = expm(A*dT);
 Bd = Ad*[0,0;
          0,0;
          1,0;
          0,1];
      
 % Specify matrices of LQ controller     
 Q = diag([1,1,0.001,0.001]);
 R = 100*diag([10000,10000]);
 
