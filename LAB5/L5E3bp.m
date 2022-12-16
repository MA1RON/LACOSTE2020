clc
clear
close all

%% Assegno i dati
%SS
rho_ss=7800; %kg/m^3
kss=15; %W/(m*K)
cp_ss=500; %J/(kg*K)
tss=1e-3; %m

%Isolante
rho_ins=2000; %kg/m^3
kins=1; %W/(m*K)
cp_ins=1000; %J/(kg*K)
tins=2e-3;  %m

Tbottom=70; %K
Tair=20; %K
hh=15; %W/(m^2*K)

%% Soluzione numerica
%Creazione della griglia
dx=1.5e-5;
xss=(0:dx:tss)';
nss=length(xss);

xins=(tss+dx:dx:tins+tss)';
xx=[xss;xins];
nn=length(xx);

%Definizione delle proprietà
kk=kss*(xx<=tss)+kins*(xx>tss);

%Assemblaggio della matrice A
sub_diag=ones(nn,1);
main=-2*ones(nn,1);
super_diag=sub_diag;

BB=[sub_diag,main,super_diag];

AA=spdiags(BB,-1:1,nn,nn);

%Assemblaggio vettore b
bb=zeros(nn,1);

%Condizioni al contorno
%i=1
AA(1,1)=1;
AA(1,2)=0;
bb(1)=Tbottom;

%i=nss
AA(nss,nss-1)=kss/kins;
AA(nss,nss)=-1-kss/kins;
AA(nss,nss+1)=1;
bb(nss)=0;

%i=end
AA(end,end-1)=-1;
AA(end,end)=1+hh*dx/kins;
bb(end)=hh*dx*Tair/kins;

%% Soluzione
TT=AA\bb;

figure(1)
plot(xx*1e3,TT,'linewidth',2)
title('Distribuzione di temperatura')
xlabel('Spessore (mm)')
ylabel('Temperaturq (°C)')
set(gca,'fontsize',18)
box on
grid on

