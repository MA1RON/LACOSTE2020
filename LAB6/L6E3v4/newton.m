function [out,res,err,jj] = newton(fun,jac,inp,tol,maxjj)
err = 1; res = 1;
out = inp; outp = out;
for jj = 1:maxjj
    out = outp - jac(outp)\fun(outp);
    res = abs(fun(out));
    err = norm(out-outp)/norm(out-inp);
    outp = out;
    
    if norm(res)<tol && err<tol
        return
    end
end
end