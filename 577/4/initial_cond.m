   function f1 = initial_cond(x)

   n=length(x);
   x=x-0.5;
   for i=1:n
     if(x(i)<0.25) f1(i)=0; end
     if(x(i)>=0.25 && x(i)<0.5) f1(i)=(x(i)-0.25)/0.25; end
     if(x(i)>=0.5  && x(i)<0.75) f1(i)=1; end
     if(x(i)>=0.75 && x(i)<1.0) f1(i)=1-(x(i)-0.75)/0.25; end
     if(x(i)>=1.0) f1(i)=0; end
   end

   f1=f1';

