clear all; close all;
%assegno i dati:
LL=1;
kk=1;cp=1;rho=1;
Tleft=0;
Tright=5;
 
%creo la griglia:
Nnodi=101;
xx=[linspace(0,LL,Nnodi)]';
deltax=LL/(Nnodi-1);
qqq=100*sin(xx*pi/LL);

TT0=linspace(Tleft,Tright,Nnodi)';
figure
plot(xx,TT0);hold on
xlabel('x(m)')
ylabel('T(^oC)')
grid on
set(gca,'Fontsize',18)

dt=0.1;
tend=2;
time=0:dt:tend;
%Assegno la matrice:
aa=ones(Nnodi,1)*dt/deltax^2*kk/cp/rho;
diag_princ=1+2*aa;
diag_sub=-aa(1:end-1);
diag_sup=diag_sub;
AA=diag(diag_princ,0)+diag(diag_sub,-1)+diag(diag_sup,1);
AA(1,1)=1;AA(1,2)=0;
AA(end,end-1)=0;AA(end,end)=1; 
 
%risolvo il sistema lineare:
for ii=2:length(time)
    bb=qqq*dt/cp/rho+TT0;
    bb(1)=Tleft;bb(end)=Tright;
    TT=AA\bb;
    TT0=TT;
    plot(xx,TT,'linewidth',3)
    pause(0.2)
end
    
 
