function fun = funbuilder(TT,kk,kder,hh,Text,qqq,dx)
nn = length(qqq);

ki = kk(TT(1));
kp = kk(TT(2:end-1));
ks = kk(TT(end));

dkp = kder(TT(2:end-1));

fun(1,1) = ki*(TT(2)-TT(1))/dx/hh-TT(1)+Text;
fun(2:nn-1,1) = dkp.*(TT(3:end)-TT(1:end-2)).^2/4/dx^2 + kp.*(TT(3:end)-2*TT(2:end-1)+TT(1:end-2))/dx^2 + qqq(2:end-1);
fun(nn,1) = ks*(TT(end-1)-TT(end))/dx/hh-TT(end)+Text;
end