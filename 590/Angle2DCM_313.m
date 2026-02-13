
function O = Angle2DCM_313(ang33,ang1,ang3)

O33 = [ cos(ang33) sin(ang33)  0;
       -sin(ang33) cos(ang33)  0;
        0          0           1];

O1 =  [ 1          0           0;
        0          cos(ang1)   sin(ang1);
        0         -sin(ang1)   cos(ang1)];

O3 =  [ cos(ang3) sin(ang3)  0;
       -sin(ang3) cos(ang3)  0;
        0         0          1];

O = O33*O1*O3;

end