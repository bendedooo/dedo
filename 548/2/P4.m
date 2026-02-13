%% PROBLEM 4
omega = 2*pi/(23*3600+56*60+4);
mu    = 398600.4405;
phi = pi/4;
r = [6378:200000];

w = (mu^2./(r.^4) + omega^4*r.^2*cos(phi)^2-2*mu*omega^2*cos(phi)^2./r);

plot(r, sqrt(w))
r(find(w==min(w)))
wmin=sqrt(min(w))
Ve=3.5;
m_0 = 100/(exp(-wmin*3600/Ve));