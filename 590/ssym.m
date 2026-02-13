function out = ssym(x)

out(1,1) = 0;
out(1,2) = -x(3);
out(1,3) = x(2);
out(2,1) = x(3);
out(2,2) = 0;
out(2,3) = -x(1);
out(3,1) = -x(2);
out(3,2) = x(1);
out(3,3) = 0;

end
