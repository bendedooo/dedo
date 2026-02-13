function [statesOut] = Adj(t,states,lambda)
tf = 10;
T = 1;

p1 = states(1);
p2 = states(2);
p1_dot = p2 * (lambda/(tf-t));
p2_dot = -p1 + p2/T;

statesOut = [p1_dot; p2_dot];
end