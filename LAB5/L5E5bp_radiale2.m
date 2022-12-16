clc
clear 
close all

%% Assegno i dati
d1=20e-2; %m
d2=23e-2; %m
LL=50e-2; %m
Tco2=300; %°C
hh=50; %W/(m^2*K)
qq=80e3; %m^3

kbarra=350; %W/(mK)
kss=50; %W/(mK)

Across1=(pi*d1^2/4); %m^2
V1=LL*Across1; %m^3

Across2=pi*(d2^2/4-d1^2/4); %m^2
V2=LL/2*Across2; %m^3

VV=V1+V2; %m^3
As=pi*d1^2/4+pi*d2^2/4; %m^2

%% Soluzione numerica
%Creazione della griglia
dr=1e-4;
r1=(0:dr:d1/2)';
n1=length(r1);

r2=(d1/2+dr:dr:d2/2)';

rr=[r1;r2];
nn=length(rr);

%Definizione delle proprietà
kk=kbarra*(rr<=d1/2)+kss*(rr>d1/2);
qvol=qq*(rr<=d1/2)+0*(rr>d1/2);

%Assemblaggio matrice A
CC=As*hh/VV./kk;
sub_diag=1-dr./(2*rr);
principale=-2-CC*dr^2;
super_diag=1+dr./(2*rr);

BB=[[sub_diag(2:end);0],principale,[0;super_diag(1:end-1)]];

AA=spdiags(BB,-1:1,nn,nn);

%Vettore dei termini noti b
bb=-qvol.*dr^2./kk-CC*dr^2*Tco2;

%Condizioni al contorno
%i=1
AA(1,1)=1;
AA(1,2)=-1;
bb(1)=0;

%i=n1
AA(n1,n1-1)=kbarra/kss;
AA(n1,n1)=-(1+kbarra/kss);
AA(n1,n1+1)=1;
bb(n1)=0;

%i=end
AA(end,end-1)=-1;
AA(end,end)=hh*dr/kss+1;
bb(end)=hh*dr*Tco2/kss;

%% Soluzione
TT=AA\bb;

figure(1)
plot(rr*1e2,TT,'linewidth',2)
grid on
box on
title ('Distribuzione di temperatura')
xlabel('Raggio (cm)')
ylabel('Temperatura (°C)')
set(gca,'fontsize',18)
