function [x0,res,err,jj] = newton(fun,der,xguess,tol,maxjj)
err = 1; res = fun(xguess); x0 = xguess; xp = x0;
for jj = 1:maxjj
   x0 = xp - fun(xp)/der(xp);
   res = abs(fun(x0));
   err = abs(x0-xp)/abs(x0-xguess);
   xp = x0;
   
   if err<tol && res<tol
       return
   end
end
end