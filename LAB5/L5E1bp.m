clc
clear all
close all

%% Assegno i dati
ti=10e-2; %m
ki=0.04; %W/mK
tb=25e-2; %m
kb=0.55; %W/mK

Tin=20; %K
Text=-5; %K
hh=25; %W/m^2K

%% Soluzione numerica
%Creazione della griglia
dxi=1e-4; 
xi=(0:dxi:ti)'; 
ni=length(xi);  %strato di isolante

dxb=1e-3;
xb=(ti+dxb:dxb:tb+ti)'; %strato di mattoni
nb=length(xb);

xx=[xi;xb]; %spessore completo
nn=length(xx); %numero di nodi totali

%Assemblaggio della matrice A
sub_diag=ones(nn,1);
main=-2*ones(nn,1);
super_diag=sub_diag;

BB=[sub_diag,main,super_diag];

AA=spdiags(BB,-1:1,nn,nn);

%Assemblaggio della matrice b
bb=zeros(nn,1);

%Condizioni al contorno
%i=1
AA(1,1)=1+hh*dxi/ki;
AA(1,2)=-1;
bb(1)=hh*dxi*Text/ki;

%i=ni condizione all'interfaccia
AA(ni,ni+1)=1;
AA(ni,ni)=-1-ki/dxi*dxb/kb;
AA(ni,ni-1)=ki/dxi*dxb/kb;
bb(ni)=0;

%i=end
AA(end,end-1)=0;
AA(end,end)=1;
bb(end)=Tin;

%Soluzione
TT=AA\bb;

figure(1)
plot(xx*1e2,TT,'linewidth',2)
hold on
plot([10 10],[min(TT) max(TT)],'--k','linewidth',2)
title('Distribuzione di temperatura')
xlabel('Spessore (cm)')
ylabel('Temperatura (°C)')
set(gca,'fontsize',22)
box on
grid on
legend('Temperatura','Interfaccia','location','se')

%Perdite termiche sulla parete
kk=linspace(ki,kb,10);

figure(3)
plot([10 10],[min(TT) max(TT)],'--k','linewidth',2)
hold on
legend('Interfaccia','location','se')

for ii=1:length(kk)
    
    AA(1,1)=1+hh*dxi/kk(ii);
    AA(1,2)=-1;
    bb(1)=hh*dxi*Text/kk(ii);
    
    AA(ni,ni+1)=1;
    AA(ni,ni)=-1-kk(ii)/dxi*dxb/kb;
    AA(ni,ni-1)=kk(ii)/dxi*dxb/kb;
    bb(ni)=0;
    
    AA(end,end-1)=0;
    AA(end,end)=1;
    bb(end)=Tin;
    
    TT=AA\bb;
    
    qlosses(ii)=hh*(TT(1)-Text);
    
    figure(3)
    plot(xx*1e2,TT,'displayname',['k = ' num2str(kk(ii)) ' W/(m*K)'],'linewidth',2)
    hold on
    legend('-dynamiclegend','location','se')

end

figure(2)
plot(kk,qlosses,'r-o','linewidth',2)
title('Variazione delle perdite sulla parete')
xlabel('Conducibilità termica (W/(m*K))')
ylabel('Perdite termiche (W/m^2)')
set(gca,'fontsize',18)
box on
grid on

figure(3)
title('Distribuzione di temperatura')
xlabel('Spessore (cm)')
ylabel('Temperatura (°C)')
set(gca,'fontsize',18)
box on
grid on

figure(4)
plot(xx*1e2,TT,'linewidth',2)
hold on
plot([10 10],[min(TT) max(TT)],'--k','linewidth',2)
title('Distribuzione di temperatura')
xlabel('Spessore (cm)')
ylabel('Temperatura (°C)')
set(gca,'fontsize',22)
box on
grid on
legend('Temperatura','Interfaccia','location','se')