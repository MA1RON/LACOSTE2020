function [yy,time]=ES05_CN(t0,t1,y0,fun,hh)
 
tic;
tt=t0:hh:t1;
nn=length(tt);
yy=ones(nn,1)*y0;
 
for jj=2:nn
    fn=fun(tt(jj-1));
    fn1=fun(tt(jj));
    yy(jj)=(hh/2*fn+hh/2*fn1+hh/2*yy(jj-1)+yy(jj-1))/(1-hh/2);
end
time=toc;
