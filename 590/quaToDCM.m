function out = quaToDCM(q)
angle = acos(q(4)).*2;
axis1 = q(1)./sin(angle/2);
axis2 = q(2)./sin(angle/2);
axis3 = q(3)./sin(angle/2);
axis = [axis1 axis2 axis3]';
% plot(angle*180/pi);
axis  = axis./norm(axis);
R = cos(angle).*eye(3) + (1-cos(angle)).*axis*axis' + sin(angle).*ssym(axis);
% orthonomalization!!!!
out = R';
end