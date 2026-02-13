function O_B_I = qToDCM(q)
% quat2rotm([q_0(4) q_0(1:3)'])

q_v = q(1:3);
q_4 = q(4);
O_B_I = (q_4^2*eye(3) - diag(sum(q_v.^2)*ones(1,3))) + 2.*q_v*q_v' - 2.*q_4*ssym(q_v); % eqn (12-13b) Wertz

% O_B_I = (q_4^2 - q_v*q_v')*eye(3) + 2.*q_v*q_v' - 2.*q_4*ssym(q_v);

end