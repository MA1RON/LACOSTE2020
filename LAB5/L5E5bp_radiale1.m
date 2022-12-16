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
As=pi*d1^2/4+pi*d2^2/4; %m^2

%% Soluzione numerica
%Creazione della griglia
dr=1e-4;
rr=(0:dr:d1/2)';
nn=length(rr);

%Assemblaggio matrice A
CC=As*hh/VV/kbarra;
sub_diag=1-dr./(2*rr);
principale=(-2-CC*dr^2)*ones(nn,1);
super_diag=1+dr./(2*rr);

BB=[[sub_diag(2:end);0],principale,[0;super_diag(1:end-1)]];

AA=spdiags(BB,-1:1,nn,nn);

%Vettore dei termini noti b
bb=(-qvol*dr^2/kbarra-CC*dr^2*Tco2)*ones(nn,1);

%Condizioni al contorno
%i=1
AA(1,1)=1;
AA(1,2)=-1;
bb(1)=0;

%i=end
AA(end,end-1)=-1;
AA(end,end)=hh*dr/kbarra+1;
bb(end)=hh*dr*Tco2/kbarra;

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
