clc
clear
close all

%% Dati
kk=0.2; %W/mK
din=20e-2; %m
dout=30e-2; %m

Taria=5; %°C
haria=100; %W/m^2K

qsurf=50; %W/m^2

qvol=0; %W/m^3

%% Griglia
dr=1e-3;
rr=(din/2:dr:dout/2)';
nn=length(rr);

%% Matrice A
sub_diag=1-dr./(rr);
main=-2*ones(nn,1);
super_diag=1+dr./(rr);

%Coordinate sferiche: manipolazione delle diagonali con gli elementi
%fittizi
BB=[[sub_diag(2:end);0],main,[0;super_diag(1:end-1)]];

AA=spdiags(BB,-1:1,nn,nn);

%% Vettore b
bb=-qvol/kk*dr^2*ones(nn,1);

%% Condizioni al contorno
%Nodo 1
AA(1,1)=1;
AA(1,2)=-1;
bb(1)=dr/kk*qsurf;

%Nodo end
AA(end,end-1)=-1;
AA(end,end)=1+haria*dr/kk;
bb(end)=haria*dr*Taria/kk;

%% Soluzione
TT=AA\bb;

%% plot
figure(1)
plot(rr*1e2,TT,'linewidth',4)
title('Distribuzione di temperatura lungo r')
xlabel('Raggio (cm)')
ylabel('Temperatura (°C)')
set(gca,'fontsize',22)
grid on
