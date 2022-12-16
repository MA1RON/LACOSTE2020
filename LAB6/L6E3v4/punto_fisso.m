function [TT,res,err,jj] = punto_fisso(fun,T0,tol,maxjj)
err = 1; res = fun(T0);
TT = T0; Tp = T0;
for jj = 1:maxjj
    TT = fun(Tp);
    res = fun(TT)-TT;
    err = norm(TT-Tp)/norm(TT-T0);
    Tp = TT;
    if err<tol && norm(res)<tol
        return
    end
end
end