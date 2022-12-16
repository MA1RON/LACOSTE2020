clc
clear
close all

%% Dati 
kk=17; %W/mK
din=25e-3; %m
e=35e-3; %m 
s=30e-3; %m

Aesagono=3*sqrt(3)*(e/2)^2/2; 
dout=sqrt(4*Aesagono/pi); %diametro equivalente: ricavato dall'area dell'esagono
                            %la sezione della corona esagonale e della
                            %corona circolare devono essere uguali affinchè
                            %la nuova geometria sia equivalente a quella di
                            %partenza

Tin=450; %K
% WTF hin vale 1000 W/m^2K
hin=500; %W/m^2K

Tout=298.15; %K
hout=100; %W/m^2K
qvol=0; %W/m^3

%% Griglia
dr=1e-6;
rr=(din/2:dr:dout/2)';
nn=length(rr);

%% Matrice A
sub_diag=1-dr./(2*rr);
main=-2*ones(nn,1);
super_diag=1+dr./(2*rr);

BB=[[sub_diag(2:end);0],main,[0;super_diag(1:end-1)]];

AA=spdiags(BB,-1:1,nn,nn);

%% Vettore b
bb=-qvol/kk*dr^2*ones(nn,1);

%% Condizioni al contorno
%Nodo 1
AA(1,1)=1+hin*dr/kk;
AA(1,2)=-1;
bb(1)=hin*Tin*dr/kk;

%Nodo end
AA(end,end-1)=-1;
AA(end,end)=1+hout*dr/kk;
bb(end)=hout*Tout*dr/kk;

%% Soluzione
TT=AA\bb;

%% plot
figure(1)
plot(rr*1e3,TT,'linewidth',4)
title('Distribuzione di temperatura lungo r')
xlabel('Raggio (mm)')
ylabel('Temperatura (K)')
set(gca,'fontsize',22)
grid on
xlim([din/2*1e3 dout/2*1e3])