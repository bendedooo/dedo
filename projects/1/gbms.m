function gbmsOut = gbms(W,mu,sigma,Y_0,dt)

W=W';
sizeW = size(W);
Y(1:sizeW(1),1) = Y_0;
for i=1:sizeW(2)-1
    b = mu * Y(:,i);
    h = sigma * Y(:,i);
    dW = W(:,i+1)-W(:,i);
    Y(:,i+1) = Y(:,i) + b.*dt + h.*dW;
end
gbmsOut=Y;

end