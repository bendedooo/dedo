function stop = HW1Q2outfun(x,optimValues,state)
global h
stop = false;
switch state
    case 'iter'
          h   = [h; x(1) x(2)]; 
    otherwise
end
end
