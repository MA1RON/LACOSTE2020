clear all;
close all;
%Dati;
Rout=33.7e-3/2;
Rin=Rout-2.9e-3;
kss=13;
%in questo specifico problema DEVO convertire le temperature in K, almeno
%per calcolare la conducibilià dell'elio.
Twater=100+273;
THemax=100+273;
THemin=20+273;
Pr=0.7;
Re=1e5;

%per calcolare il coefficiente di scambio termico, mi devo ricordare
% che la correlazione di Dittus Boelter richiede sempre le prop del fluido
%alla temperatura media tra ingresso e uscita: qui le temperature sono
%entrambe assgenate, quindi il problema è FINTAMENTE non lineare!
pp=5;
Tave=0.5*(THemax+THemin);
khe=@(temp) (2.682e-3*(1+1.123e-3*pp)*temp.^(0.71*(1-2e-4*pp)));
hh=0.023*Re^0.8*Pr^0.7*khe(Tave)/2/Rout;

% a questo punto il problema si risolve in geometria radiale, senza
% intoppi:
dr=1e-4;
rr=[Rin:dr:Rout]';
enne=length(rr);
main=-2*ones(enne,1);
sup=[0;ones(enne-1,1)+dr/2./rr(1:end-1)];
sub=[ones(enne-1,1)-dr/2./rr(2:end);0];
banda=[sub,main,sup];
AA=spdiags(banda,-1:1,enne,enne);
AA(end,end)=1;
AA(end,end-1)=0;
%min THe
AA(1,1)=-1-hh*dr/kss;
AA(1,2)=1;
bb=zeros(enne,1);
bb(end)=Twater;
bb(1)=-hh*dr/kss*THemin;%all'ingresso del tubo la temperatura dell'He è la minima

%risolvo il problema
TT=AA\bb;
%faccio il grafico della soluzione:
plot(rr,TT)
xlabel('Posizione radiale (mm)')
ylabel('Temperatura (K)')
grid on
