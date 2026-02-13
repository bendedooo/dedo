function accMeasOut = accMeas(g,phi,t)
accMeasOut = [ -1-g*sin(phi)*sin(t);
               -g*sin(phi)*cos(t);
               -g*cos(phi)];
end