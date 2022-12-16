clc
clear
close all

%Dati
kb=0.55; %W/mK
ki=0.04; %W/mK

tb=25e-2; %m
ti=10e-2; %m
thick=tb+ti; %m

Tin=20; %°C
Tout=-5; % °C
hout=25; %W/m^2K

kk=(kb*tb+ki*ti)/thick; %  k equivalente: media pesata sugli spessori

%% Griglia
dx=1e-4;
xx=(0:dx:thick)';
nn=length(xx);

%% Matrice A
sub_diag=ones(nn,1);
main=-2*ones(nn,1);
super_diag=sub_diag;

BB=[sub_diag,main,super_diag];

AA=spdiags(BB,-1:1,nn,nn);

%% Vettore b
bb=zeros*ones(nn,1);

%% Condizioni al contorno
%Nodo 1
AA(1,1)=1+hout*dx/kk;
AA(1,2)=-1;
bb(1)=hout*dx*Tout/kk;

%Nodo end
AA(end,end-1)=0;
AA(end,end)=1;
bb(end)=Tin;

%% Soluzione
TT=AA\bb;

%% Perdita sulla parete dell'isolante
qlosses=hout*(TT(1)-Tout); %(W/m^2)

figure(1)
plot(xx*1e2,TT,'linewidth',4)
title('Distribuzione di temperatura')
xlabel('Spessore (cm)')
ylabel('Temperatura (°C)')
set(gca,'fontsize',22)
xlim([0 thick*1e2])
grid on