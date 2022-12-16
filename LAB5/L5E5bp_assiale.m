clc
clear 
close all

%% Assegno i dati
d1=20e-2; %m
d2=23e-2; %m
LL=50e-2; %m
Tco2=300; %°C
hh=50; %W/(m^2*K)
qvol=80e3; %m^3

kbarra=350; %W/(mK)
kss=50; %W/(mK)

Across1=(pi*d1^2/4); %m^2
V1=LL*Across1; %m^3

Across2=pi*(d2^2/4-d1^2/4); %m^2
V2=LL/2*Across2; %m^3

VV=V1+V2; %m^3
As=pi*d1*LL/2+pi*d2*LL/2; %m^2

%Conducilità equivalente barra-acciaio
keq=(kbarra*Across1+kss*Across2)/(Across1+Across2);

%% Soluzione numerica
%Creazione della griglia
dx=1e-3;
x1=(0:dx:LL/2)';
n1=length(x1);

x2=(LL/2+dx:dx:LL)';
xx=[x1;x2];
nn=length(xx);

%Definizione delle proprietà
kk=kbarra*(xx<=LL/2)+keq*(xx>LL/2);

%Assemblaggio matrice A
CC=As*hh/VV./kk;
sub_diag=ones(nn,1);
principale=-2-CC*dx^2;
super_diag=ones(nn,1);

BB=[sub_diag,principale,super_diag];

AA=spdiags(BB,-1:1,nn,nn);

%Vettore dei termini noti b
bb=-qvol*dx^2./kk-CC*dx^2*Tco2;

%Condizioni al contorno
%i=1
AA(1,1)=hh*dx/kbarra+1;
AA(1,2)=-1;
bb(1)=hh*dx*Tco2/kbarra;

%i=n1
AA(n1,n1-1)=kbarra/keq;
AA(n1,n1)=-(1+kbarra/keq);
AA(n1,n1+1)=1;
bb(n1)=0;

%i=end
AA(end,end-1)=-1;
AA(end,end)=hh*dx/keq+1;
bb(end)=hh*dx*Tco2/keq;

%% Soluzione
TT=AA\bb;

figure(1)
plot(xx*1e2,TT,'linewidth',2)
grid on
box on
title ('Distribuzione di temperatura')
xlabel('Lunghezza (cm)')
ylabel('Temperatura (°C)')
set(gca,'fontsize',18)
