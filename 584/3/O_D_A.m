function O_mat = O_D_A (a,b,c)
%  orientation matrix for a O1(a)*O2(b)*O3(c) rotation
% a,b,c are either scalars or (1,n) vectors.
%returns a (3,3*n) matrix.

Cb=  cos(b);
Sb=  sin(b);
Ca=  cos(a);
Sa=  sin(a);
Cc=  cos(c);
Sc=  sin(c);
N_s = size(a,2);
O_temp =[Cb.*Cc         Cb.*Sc              -Sb;
    Cc.*Sa.*Sb-Ca.*Sc   Ca.*Cc+Sa.*Sb.*Sc   Cb.*Sa;
    Sa.*Sc+Ca.*Cc.*Sb   Ca.*Sa.*Sc-Cc.*Sa   Ca.*Cb];
O_mat = [];
for cur_i = 1: N_s
    O_mat  = [O_mat O_temp(:,cur_i:N_s:end)];
end

% O_temp(cur_i+N_s:N_s:end-N_s,:) ;O_temp(cur_i+2*N_s:N_s:end,:)
end