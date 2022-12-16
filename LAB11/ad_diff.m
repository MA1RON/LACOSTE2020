%Problema: un tubo di 1 cm di diametro interno e 1 mm di spessore è percorso 
%da una portata di 0.12 kg/s di acqua a pressione atmosferica in deflusso
%turbolento. Il condotto in acciaio, ben coibentato al suo esterno, viene
%sottoposto, da t =0; a un riscaldamento induttivo per tutta la sua lunghezza
% L = 2m. La potenza depositata per unità di lunghezza del condotto è pari a:
% Q=q0*exp(-5*x^2), dove x è la coordinata spaziale rispetto al centro del
% condotto e q0=2e3 W/m. La temperatura di ingresso dell'acqua è 25 oC.
%Utilizzando un modello 1D, calcolare il profilo spaziale della 
%temperatura nel condotto e della temperatura di bulk nel fluido a t = 1s,
%5s, 10s e in condizioni stazionarie, supponendo che il coefficiente di scambio
%termico all'interno del condotto sia costante e pari a 2000 W/m^2K.

% clear all;close all
%dati del problema:
LL=2; %[m];
q0=2e3; %[W/m];
qq=@(xx) (q0*exp(-5*(xx-1).^2)); %traslato per avere x=0 all'ingresso del condotto
portata=0.12; %[kg/s]
diamint=0.01; %[m]
spess=1e-3; %[m]
dens=1000; %[kg/m3]
cpacqua=4186; %[J/kgK]
rhoss=7800; %[kg/m3]
cpss=500; %[J/kgK]
kappass=30; %[W/mK]
acca=2000; %[W/m2K]
Tin=25; %[K]

%grandezze derivate dai dati del problema:
areaflow=pi*diamint^2/4;
speed=portata/areaflow/dens;
%gli istanti di tempo a cui mi è richiesto l'output:
timeout1=1;timeout2=5;timeout3=10;

%prima provo ad applicare la sorgente di calore al fluido direttamente
%rho*cp*(dT/dt+u*dT/dx)=qvol
%dT/dt+u*dT/dx=qvol/(rho*cp)
%schema numerico: uso un upwind con Eulero implicito
%Ti,new-Ti,old+u*deltat/deltax*(Ti,new-Ti-1,new)=deltat*qvoli/(rho*cp)
%aaa=u*deltat/deltax;
%-aaa*Ti-1,new+(1+aaa)*Ti,new=Ti,old+deltat*qvol/(rho*cp)

%parametri numerici
deltax=0.01;
nnodi=LL/deltax+1;
xx=linspace(0,LL,nnodi)';
deltat=1e-1;

%definisco i coefficienti che uso per gli elementi di matrice
aaa=speed*deltat/deltax;
%Costruisco la diagonale principale:
maindiag=(1+aaa)*ones(size(xx));
% e la SOTTOdiagonale (SOTTO e non SOPRA per via dell'upwind)
subdiag=(-aaa)*ones(size(xx));
%assemblo la matrice
AAA=spdiags([subdiag,maindiag],[-1 0],nnodi,nnodi);
AAA(1,1)=1;
qvol=qq(xx)/areaflow;
bb=deltat*qvol/dens/cpacqua;
bb(1)=0;
Told=Tin*ones(nnodi,1);

toll=1e-2;
variaz=10*toll;
time=0;ii=0;
figure
set(gca,'fontsize',24)

iplot=0;
while variaz>toll
    ii=ii+1;
    time=time+deltat;
    Tnew=AAA\(Told+bb);
    variaz=abs(Tnew(end)-Told(end));
    Told=Tnew;
    if abs(time-timeout1)<= deltat/2
        iplot=1;
    elseif abs(time-timeout2)<= deltat/2
        iplot=1;
    elseif abs(time-timeout3)<= deltat/2
        iplot=1;
    end
    if iplot==1
        plot(xx,Tnew,'linewidth',2)
        hold on
        iplot=0;
    end
    %pause(0.2)
end
plot(xx,Tnew,'linewidth',2)
legend('Time = 1 s','Steady state')
xlabel('x (m)')
ylabel('Temperature (^oC)')





%ora risolviamo il problema per il solido, considerando il fluido sempre a Tin:
%rhoss*cpss*(dT/dt)=kss*(d2T/dt2)+qvol-acca*(pi*Din*LL)/(areatube*LL)*(T-Twater)
%(dT/dt)=kss/(rhoss*cpss)*(d2T/dt2)+qvol/(rhoss*cpss)-acca*(pi*Din*LL)/(areatube*LL)/(rhoss*cpss)*(T-Twater)
%alpha=kappass/rhoss/cpss*deltat/deltax^2;
%convez=acca*(pi*Din*LL)/(areatube*LL)/(rhoss*cpss)
%Ti,new-Ti,old=alpha*(Ti-1,new-2*Ti,new+Ti+1,new)+qvol/(rhoss*cpss)*deltat-convez*deltat*(Ti,new-Tin);
%-alpha*(Ti-1,new)+(1+2*alpha+convez*deltat)*Ti,new-alpha*(Ti+1,new)=Ti,old+qvol/(rhoss*cpss)*deltat+convez*deltat*Tin;

alpha=kappass/rhoss/cpss*deltat/deltax^2;
areatubo=pi*((diamint/2+spess)^2-diamint^2/4);
convez=acca*(pi*diamint)/(areatubo)/(rhoss*cpss);
convez_acqua=acca*(pi*diamint)/(areaflow)/(dens*cpacqua);
qvolss=qq(xx)/areatubo;
maindiagss=(1+2*alpha+convez*deltat)*ones(size(xx));
subdiagss=(-alpha)*ones(size(xx));
supdiagss=subdiagss;
AAAss=spdiags([subdiagss,maindiagss,supdiagss],[-1 0 1],nnodi,nnodi);
AAAss(1,1)=1;AAAss(1,2)=-1;
AAAss(end,end)=1;AAAss(end,end-1)=-1;
bbss=(qvolss/rhoss/cpss+convez*Tin)*deltat;

Toldss=Tin*ones(nnodi,1);
toll=1e-2;
variaz1=10*toll;
time1=0;jj=0;
figure
set(gca,'fontsize',24)

iplot=0;
while variaz1>toll
    jj=jj+1;
    time1=time1+deltat;
    bbnew=Toldss+bbss;
    bbnew(1)=0;
    bbnew(end)=0;
    Tnewss=AAAss\bbnew;
    variaz1=abs(max(Tnewss)-max(Toldss));
    Toldss=Tnewss;
    if abs(time1-timeout1)<= deltat/2
        iplot=1;
    elseif abs(time1-timeout2)<= deltat/2
        iplot=1;
    elseif abs(time1-timeout3)<= deltat/2
        iplot=1;
    end
    if iplot==1
        plot(xx,Tnewss,'linewidth',2)
        hold on
        iplot=0;
    end
end
plot(xx,Tnewss,'linewidth',2)
legend('Time = 1 s','Time = 5 s','Time = 10 s','Steady state')
xlabel('x (m)')
ylabel('Temperature (^oC)')

%ora cerchiamo di mettere insime i due studi:
% faccio una griglia che abbia un nodo sul solido e il corrispondente nodo
% sul fluido:
xx1=ones(2*nnodi,1);
xx1(1:2:end,1)=xx;
xx1(2:2:end,1)=xx;

maindiagall(1:2:2*nnodi,1)=maindiagss;
maindiagall(2:2:2*nnodi,1)=maindiag+convez_acqua*deltat;
 supdiagall(2:2:2*nnodi,1)=-convez*deltat*ones(nnodi,1); %okkio al conto di spdiags
 supdiagall(1:2:2*nnodi,1)=zeros(nnodi,1);
subdiagall(2:2:2*nnodi,1)=zeros(nnodi,1);
subdiagall(1:2:2*nnodi,1)=-convez_acqua*deltat*ones(nnodi,1); %okkio al conto di spdiags
subsubdiagall(1:2:2*nnodi,1)=subdiagss;
subsubdiagall(2:2:2*nnodi,1)=subdiag;
 supsupdiagall(1:2:2*nnodi,1)=subdiagss;
 supsupdiagall(2:2:2*nnodi,1)=zeros(nnodi,1);

AAall=spdiags([subsubdiagall, subdiagall,maindiagall, supdiagall ,supsupdiagall],[-2:2],2*nnodi,2*nnodi);
bball(1:2:2*nnodi,1)=(qvolss/rhoss/cpss)*deltat;
bball(2:2:2*nnodi,1)=zeros(nnodi,1);
%%BC
AAall(1,1)=kappass;
AAall(1,3)=-kappass;
AAall(1,2)=0;
AAall(2,1)=0;
AAall(2,2)= 1;
AAall(end-1,end-1)=kappass;
AAall(end-1,end)=0;
AAall(end-1,end-3)=-kappass;
AAall(end,end)=1+aaa;
AAall(end,end-1)=0;
AAall(end,end-2)=-aaa;

%%soluzione


Toldall=Tin*ones(2*nnodi,1);
toll=1e-2;
variaz2=10*toll;
time2=0;zz=0;
figure
set(gca,'fontsize',24)
iplot=0;
while variaz2>toll
    zz=zz+1;
    time2=time2+deltat;
    bbnewall=Toldall+bball;
    bbnewall(1,1)=0;bbnew(2,1)=Tin;
    bbnewall(end-1,1)=0;
    Tnewall=AAall\bbnewall;
    Tsolid=Tnewall(1:2:end,1);
    Tfluid=Tnewall(2:2:end,1);
    variaz2=max(abs(max(Tsolid)-max(Toldall(1:2:end))),(Tfluid(end)-Toldall(end)));
    Toldall=Tnewall;
    if abs(time2-timeout1)<= deltat/2
        iplot=1;
    elseif abs(time2-timeout2)<= deltat/2
        iplot=1;
    elseif abs(time2-timeout3)<= deltat/2
        iplot=1;
    end
    if iplot==1
        plot(xx,Tsolid,'--','linewidth',2)
        hold on
        plot(xx,Tfluid,'-','linewidth',2)
        iplot=0;
    end
end
plot(xx,Tsolid,'--','linewidth',2)
 hold on
 plot(xx,Tfluid,'-','linewidth',2)
legend('Time = 1 s, SS','Time = 1 s, H20',...
    'Time = 5 s, SS','Time = 5 s, H20',...
    'Time = 10 s, SS','Time = 10 s, H20','Steady state SS','Steady state H20')
xlabel('x (m)')
ylabel('Temperature (^oC)')



