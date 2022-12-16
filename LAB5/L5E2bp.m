clc
clear 
close all

%% Assegno i dati
%Cu
kcu=150; %W/mK
cpcu=350; %J/kgK
rhocu=8900; %kg/m^3
rhoel_cu=2e-8; %ohm*m

%Ins
kins=0.5; %W/mK
cpins=1800; %J/kgK
rhoins=2500; %kg/m^3

%N2
Tn2=77; %K
hint=1200; %W/m^2K

%Air
Ta=25+273.15; %K 
hext=15; %W/m^2K

II=20e3; %A

d1=30e-3; %m
d2=50e-3; %m
d3=52e-3; %m

L=1; %m 
S=pi*(d2^2-d1^2)/4; %m^2
V=S*L; %m^3

qcu=rhoel_cu*II^2*L/S/V; %W/m^3

%% Soluzione numerica (1)
%Creazione della griglia
drcu=1e-4;
rcu=(d1/2:drcu:d2/2)'; 
ncu=length(rcu); %strato di rame

drins=1e-5;
rins=(d2/2+drins:drins:d3/2)';
nins=length(rins); %strato di isolante

rr=[rcu;rins]; %geometria intera
nn=length(rr); %numero di nodi totali

%Definizione delle proprietà
kk=kcu*(rr<=d2/2)+kins*(rr>d2/2); %vettore delle conducibilità
qvol=qcu*(rr<=d2/2)+0*(rr>d2/2); %vettore del qvol
dr=drcu*(rr<=d2/2)+drins*(rr>d2/2); %vettore dei dr

%Assemblaggio matrice A
sub_diag=1-dr./(2*rr);
main=-2*ones(nn,1);
super_diag=1+dr./(2*rr);

BB=[[sub_diag(2:end);0],main,[0;super_diag(1:end-1)]];

AA=spdiags(BB,-1:1,nn,nn);

%Assemblaggio vettore b
bb=-qvol.*dr.^2./kk;

%Condizioni al contorno
%i=1
AA(1,1)=1+hint*drcu/kcu;
AA(1,2)=-1;
bb(1)=hint*drcu*Tn2/kcu;

%i=ncu
AA(ncu,ncu-1)=kcu/kins*drins/drcu;
AA(ncu,ncu)=-1-kcu/kins*drins/drcu;
AA(ncu,ncu+1)=1;
bb(ncu)=0;

%i=end
AA(end,end-1)=-1;
AA(end,end)=1+hext*drins/kins;
bb(end)=hext*drins*Ta/kins;

%Soluzione
TT=AA\bb;

%Grafico
figure(1)
plot(rr*1e3,TT,'linewidth',3)
hold on
plot([25 25],[135 144],'--k','linewidth',3)
title('Distribuzione di temperatura')
xlabel('Raggio (mm)')
ylabel('Temperatura (K)')
set(gca,'fontsize',22)
box on
grid on
legend('Temperatura','Interfaccia','location','nw')
xlim([min(rr*1e3) max(rr*1e3)])
ylim([135 144])

%% Soluzione numerica (2) - costo computazionale
NP=[1e2 2e2 5e2 1e3 2e3  5e3 1e4]; %numero di nodi
tt=size(rr);

for ii=1:length(NP)

    drcu = (d2/2-d1/2)/NP(ii); % passo rame
    rcu=(d1/2:drcu:d2/2)';
    ncu=length(rcu);
    
    drins = (d3/2-d2/2)/NP(ii); % passo isolante
    rins=(d2/2+drins:drins:d3/2)';
    rr=[rcu;rins];
    nn(ii)=length(rr);

    kk=kcu*(rr<=d2/2)+kins*(rr>d2/2);
    qvol=qcu*(rr<=d2/2)+0*(rr>d2/2);
    dr=drcu*(rr<=d2/2)+drins*(rr>d2/2);
    
    sub_diag=1-dr./(2*rr);
    main=-2*ones(nn(ii),1);
    super_diag=1+dr./(2*rr);
    
    BB=[[sub_diag(2:end);0],main,[0;super_diag(1:end-1)]];
    
    AA=spdiags(BB,-1:1,nn(ii),nn(ii));
    
    bb=-qvol.*dr.^2./kk;
    
    AA(1,1)=1+hint*drcu/kcu;
    AA(1,2)=-1;
    bb(1)=hint*drcu*Tn2/kcu;
    
    AA(ncu,ncu-1)=kcu/kins*drins/drcu;
    AA(ncu,ncu)=-1-kcu/kins*drins/drcu;
    AA(ncu,ncu+1)=1;
    bb(ncu)=0;
    
    AA(end,end-1)=-1;
    AA(end,end)=1+hext*drins/kins;
    bb(end)=hext*drins*Ta/kins;
    
    tic %tempo di inizio
    
    TT=AA\bb;
    
    tt(ii)=toc; %tempo di fine
 
end

figure(2)
semilogy(nn,tt*1e3,'-o','linewidth',3)
title('Costo computazionale in funzione delle dimensioni della matrice')
xlabel('#Nodi (-)')
ylabel('Costo computazionale (ms)')
set(gca,'fontsize',22)
box on
grid on

%% Soluzione numerica (3) - Raggio critico

tc=linspace(0.1e-3,15e-3,30);
drcu=1e-6;
drins=1e-7;
qlosses=size(tc);

for jj=1:length(tc)
    
    rin=d2/2+tc(jj);
    d3=2*rin;
    
    rcu=(d1/2:drcu:d2/2)';
    ncu=length(rcu);
    
    rins=(d2/2+drins:drins:d3/2)';
    rr=[rcu;rins];
    nn=length(rr);
    
    kk=kcu*(rr<=d2/2)+kins*(rr>d2/2);
    qvol=qcu*(rr<=d2/2)+0*(rr>d2/2);
    dr=drcu*(rr<=d2/2)+drins*(rr>d2/2);
    
    sub_diag=1-dr./(2*rr);
    main=-2*ones(nn,1);
    super_diag=1+dr./(2*rr);
    
    BB=[[sub_diag(2:end);0],main,[0;super_diag(1:end-1)]];
    
    AA=spdiags(BB,-1:1,nn,nn);
    
    bb=-qvol.*dr.^2./kk;
    
    AA(1,1)=1+hint*drcu/kcu;
    AA(1,2)=-1;
    bb(1)=hint*drcu*Tn2/kcu;
    
    AA(ncu,ncu-1)=kcu/kins;
    AA(ncu,ncu)=-1-kcu/kins;
    AA(ncu,ncu+1)=1;
    bb(ncu)=0;
    
    AA(end,end-1)=-1;
    AA(end,end)=1+hext*drins/kins;
    bb(end)=hext*drins*Ta/kins;
    
    TT=AA\bb;
    
    qlosses(jj)=abs(pi*d3*hext*(TT(end)-Ta)); %perdite termiche

end

[Qmax,imax]=max(qlosses);
tcrit=tc(imax);
T_theo=kins/hext-d2/2;
qloss_t=abs(pi*(d2+2*T_theo)*hext*(TT(end)-Ta));

figure(3)
plot(tc*1e3,qlosses,'r-o','linewidth',3)
hold on
plot(tcrit*1e3,Qmax,'ko','markerfacecolor','y','markersize',13)
plot([T_theo*1e3 T_theo*1e3],[min(qlosses) max(qlosses)+2],'k--','linewidth',3)
title('Perdite termiche attraverso il materiale isolante')
xlabel('Spessore_{ins} (mm)')
ylabel('Perdite termiche (W/m)')
set(gca,'fontsize',22)
box on
grid on
legend('Perdite termiche','Raggio critico numerico',...
    'Raggio critico teorico','location','se')
ylim([min(qlosses) max(qlosses)+2])