function out = ssym4(x)

% for [qv q]
out(1,1) = 0;
out(1,2) = x(3);
out(1,3) = -x(2);
out(1,4) = x(1);

out(2,1) = -x(3);
out(2,2) = 0;
out(2,3) = x(1);
out(2,4) = x(2);


out(3,1) = x(2);
out(3,2) = -x(1);
out(3,3) = 0;
out(3,4) = x(3);

out(4,1) = -x(1);
out(4,2) = -x(2);
out(4,3) = -x(3);
out(4,4) = 0;


% for [q qv]
% out(1,1) = 0;
% out(1,2) = -x(1);
% out(1,3) = -x(2);
% out(1,4) = -x(3);
% 
% out(2,1) = x(1);
% out(2,2) = 0;
% out(2,3) = x(3);
% out(2,4) = -x(2);
% 
% 
% out(3,1) = x(2);
% out(3,2) = -x(3);
% out(3,3) = 0;
% out(3,4) = x(1);
% 
% out(4,1) = x(3);
% out(4,2) = x(2);
% out(4,3) = -x(1);
% out(4,4) = 0;

end