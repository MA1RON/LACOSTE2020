function fun = funbuilder_fp(TT,kk,kder,hh,Text,qqq,dx)
nn = length(qqq);

ki = kk(TT(1)); aa = ki/hh/dx;
kp = kk(TT(2:end-1));
ks = kk(TT(end));

dkp = kder(TT(2:end-1));

fun(1,1) = (aa*TT(2)+Text)/(1+aa);
fun(2:nn-1,1) = (dkp.*(TT(3:end)-TT(1:end-2)).^2/4./kp + (TT(3:end)+TT(1:end-2)) + qqq(2:end-1)*dx^2./kp)/2;
fun(nn,1) = (aa*TT(end-1)+Text)/(1+aa);
end