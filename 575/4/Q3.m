syms t c tp mt ub t gm
expr = c/(-t+tp+mt*c/ub);
F = int(expr,t)
int(expr,t,[0 tp])


expr2 = c*log(1+t/(mt*c/ub)) -gm*t
F2 = int(expr2,t)