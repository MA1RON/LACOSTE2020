clc
clear 
close all

%% Assegno i dati 
r1=10e-3; %m
r2=14e-3; %m
r3=17e-3; %m
r4=18.5e-3; %m

Tin=25; %K
hin=2e3; %W/(m^2*K)
 
Tout=25; %K
hout=5; %W/(m^2*K)

Inom=16e3; %A

%SS
kss=1.5; %W/(m*K)
rhoelss=5e-6; %ohm*m
Rss=rhoelss/(pi*(r2^2-r1^2)); 

%Cu
kcu=30; %W/(m*K)
rhoelcu=1.8e-8; %ohm*m
Rcu=rhoelcu/(pi*(r3^2-r2^2));

%Ins
kins=0.21; %W/(m*K)
rhoelins=1e-2; %ohm*m
Rins=rhoelins/(pi*(r4^2-r3^2));

%Partitore di corrente
%Corrente nell'acciaio
Iss=Inom*(Rcu*Rins)/(Rss*Rcu+Rcu*Rins+Rins*Rss); %A
qss=Iss^2*rhoelss/(pi*(r2^2-r1^2))^2; %W/m^3

%Corrente nel rame
Icu=Inom*(Rins*Rss)/(Rss*Rcu+Rcu*Rins+Rins*Rss); %A
qcu=Icu^2*rhoelcu/(pi*(r3^2-r2^2))^2; %W/m^3

%Corrente nell'isolante
Iins=Inom*(Rcu*Rss)/(Rss*Rcu+Rcu*Rins+Rins*Rss); %A
qins=Iins^2*rhoelins/(pi*(r4^2-r3^2))^2; %W/m^3

%% Soluzione numerica
%Creazione della griglia
dr=1e-4; 
rss=(r1:dr:r2)';
nss=length(rss);

rcu=(r2+dr:dr:r3)';
ncu=nss+length(rcu);

rins=(r3+dr:dr:r4)';

rr=[rss;rcu;rins];
nn=length(rr);

%Definizione delle proprietà
kk=kss*(rr<=r2)+kcu*(rr>r2 & rr<=r3)+kins*(rr>r3);
qnom=qss*(rr<=r2)+qcu*(rr>r2 & rr<=r3)+qins*(rr>r3);

%Assemblaggio della matrice A
sub_diag=1-dr./(2*rr);
main=-2*ones(nn,1);
super_diag=1+dr./(2*rr);

BB=[[sub_diag(2:end);0],main,[0;super_diag(1:end-1)]];

AA=spdiags(BB,-1:1,nn,nn);

%Assemblaggio del vettore b
bb=-qnom.*dr.^2./kk;

%Condizioni al contorno
%i=1
AA(1,1)=1+hin*dr/kss;
AA(1,2)=-1;
bb(1)=hin*dr/kss*Tin;

%i=nss
AA(nss,nss-1)=kss/kcu;
AA(nss,nss)=-(1+kss/kcu);
AA(nss,nss+1)=1;
bb(nss)=0;

%i=ncu
AA(ncu,ncu-1)=kcu/kins;
AA(ncu,ncu)=-(1+kcu/kins);
AA(ncu,ncu+1)=1;
bb(ncu)=0;

%i=end
AA(end,end-1)=-1;
AA(end,end)=1+hout*dr/kins;
bb(end)=hout*dr/kins*Tout;

%% Soluzione
TT=AA\bb;

figure(1)
plot(rr*1e3,TT,'linewidth',2)
grid on
box on
title('Distribuzione di temperatura')
xlabel('Raggio (mm)')
ylabel('Temperatura (°C)')
set(gca,'fontsize',22)