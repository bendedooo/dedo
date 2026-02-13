function O = Angle2DCM_321(phi,theta,psi)

O1 = [1           0            0;
      0    cos(phi)     sin(phi);
      0   -sin(phi)    cos(phi)];

O2 = [cos(theta)  0  -sin(theta);
      0           1            0; 
      sin(theta)  0  cos(theta)];

O3 = [ cos(psi)   sin(psi)     0; 
      -sin(psi)   cos(psi)     0;
              0          0     1];

O = O1*O2*O3;

end