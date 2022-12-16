clc
clear 
close all

%% Assegno i dati
LL=50e-2; %m
bb=10e-2; %m
th=2e-2; %m
kk=350; %W/(m*K)
hh=500; %W/(m^2*K)
The=4.5; %K
II=7e3; %A
rho_el0=1.75e-8; %ohm*m
L0=56e-2; %m

rhoel=@(x) rho_el0*(cos(pi*x./L0)).^2; %ohm*m
As=2*(bb+th)*LL; %m^2
SS=th*bb; %m^2
VV=SS*LL; %m^3

%% Soluzione numerica
%Creazione della griglia
dx=1e-4;
xx=(0:dx:LL/2)';
nn=length(xx);

%Assemblaggio matrice A
CC=hh*As/VV/kk;

sub_diag=ones(nn,1);
main=-(2+CC*dx^2)*ones(nn,1);
super_diag=sub_diag;

BB=[sub_diag,main,super_diag];

AA=spdiags(BB,-1:1,nn,nn);

%Vettore dei termini noti b
qvol=II^2*LL*rhoel(xx)./(SS*VV); %W/m^3
bb=-qvol*dx^2./kk-CC*dx^2*The;

%Condizioni al contorno
%i=1
AA(1,1)=1;
AA(1,2)=-1;
bb(1)=0;

%i=end
AA(end,end-1)=0;
AA(end,end)=1;
bb(end)=The;

%% Soluzione
TT=AA\bb;

figure(1)
plot(xx,TT,'linewidth',2)
title('Distribuzione di temperatura')
xlabel('Lunghezza (m)')
ylabel('Temperatura (K)')
set(gca,'fontsize',18)
box on
grid on
