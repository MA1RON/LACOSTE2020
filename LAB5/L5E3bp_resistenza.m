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

%Guarnizione, ipotizzo che la resistenza di contatto possa essere
%approssimata al rapporto tra lo spessore e la condubilità, così ricavo la
%conducibilità della guarnizione
Rguarn=2.5e-4; %(K*m^2)/W
tguarn=0.01e-3; %m
kguarn=tguarn/Rguarn; %W/(m*K)

Tbottom=70; %K
Tair=20; %K
hh=15; %W/(m^2*K)

%% Soluzione numerica
%Creazione della griglia
dx=1e-7;
xss=(0:dx:tss)';
nss=length(xss);

xguarn=(tss+dx:dx:tss+tguarn)';
nguarn=nss+length(xguarn);

xins=(tss+tguarn+dx:dx:tins+tss+tguarn)';
xx=[xss;xguarn;xins];
nn=length(xx);

%Definizione delle proprietà
kk=kss*(xx<=tss)+kguarn*(xx>tss & xx<=(tss+tguarn))+kins*(xx>(tss+tguarn));

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
AA(nss,nss-1)=kss/kguarn;
AA(nss,nss)=-1-kss/kguarn;
AA(nss,nss+1)=1;
bb(nss)=0;

%i=nguarn
AA(nguarn,nguarn-1)=kguarn/kins;
AA(nguarn,nguarn)=-1-kguarn/kins;
AA(nguarn,nguarn+1)=1;
bb(nguarn)=0;

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
ylabel('Temperatura (°C)')
set(gca,'fontsize',18)
box on
grid on
xlim([xx(1)*1e3 xx(end)*1e3])

%Conservazione dell'energia
As=0.6*0.6; %m^2
Rss=tss/kss/As; %K/W
Rg=Rguarn/As; %K/W
Rins=tins/kins/As; %K/W
R=Rss+Rins+Rg; %K/W
U=1/R/As; %W/(m^2*K)

q1=abs(hh*(TT(end)-Tair));
q2=abs(U*(Tbottom-TT(end)));

