function [out,res,err,jj] = bisezione(fun,xguess,tol,maxjj)
err = 1; res = fun(xguess); x0 = xguess; xp = x0;
for jj = 1:maxjj
   if fun((x0(1)+x0(2))/2)*fun(x0(1)) <= 0
       x0(2) = (x0(1)+x0(2))/2;
   elseif fun((x0(1)+x0(2))/2)*fun(x0(2)) <= 0
       x0(1) = (x0(1)+x0(2))/2;
   else
       fprintf('bisezione non converge.\n')
       return
   end
   res = norm(fun(x0));
   err = norm(x0-xp)/norm(x0-xguess);
   xp = x0;
   
   if err<tol && res<tol
       return
   end
   out = ((x0(1)+x0(2))/2);
end
end