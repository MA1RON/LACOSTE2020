function [yy,time]=ES05_FE(t0,t1,y0,fun,hh)
 
tic;
tt=t0:hh:t1;
nn=length(tt);
yy=ones(nn,1)*y0;
for jj=2:nn
    fn=fun(tt(jj-1))+yy(jj-1);
    yy(jj)=hh*fn+yy(jj-1);
end
time=toc;
