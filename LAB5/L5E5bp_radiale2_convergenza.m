clc
clear 
close all

%% Assegno i dati
d1=20e-2; %m
d2=23e-2; %m
LL=50e-2; %m
Tco2=300; %°C
hh=50; %W/(m^2*K)
qq=80e3; %m^3

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
drv=[1e-5 5e-5 1e-4 5e-4 1e-3];
err=size(drv);
RA=err;

for ii=1:length(drv)
    
    dr=drv(ii);
    r1=(0:dr:d1/2)';
    n1=length(r1);
    
    r2=(d1/2+dr:dr:d2/2)';
    
    rr=[r1;r2];
    nn=length(rr);
    
    kk=kbarra*(rr<=d1/2)+kss*(rr>d1/2);
    qvol=qq*(rr<=d1/2)+0*(rr>d1/2);
    
    CC=As*hh/VV./kk;
    sub_diag=1-dr./(2*rr);
    principale=-2-CC*dr^2;
    super_diag=1+dr./(2*rr);
    
    BB=[[sub_diag(2:end);0],principale,[0;super_diag(1:end-1)]];
    
    AA=spdiags(BB,-1:1,nn,nn);
    
    bb=-qvol.*dr^2./kk-CC*dr^2*Tco2;
    
    AA(1,1)=1;
    AA(1,2)=-1;
    bb(1)=0;
    
    AA(n1,n1-1)=kbarra/kss;
    AA(n1,n1)=-(1+kbarra/kss);
    AA(n1,n1+1)=1;
    bb(n1)=0;
    
    AA(end,end-1)=-1;
    AA(end,end)=hh*dr/kss+1;
    bb(end)=hh*dr*Tco2/kss;
    
    TT=AA\bb;

    figure(1)
    plot(rr,TT,'displayname',['dr=' num2str(drv(ii)) ' m'],...
        'linewidth',2)
    hold on
    legend('-dynamiclegend','location','sw')
    
    KK=condest(AA);
    
    if ii==1
        drref=dr;
        Tref=TT;
        
        err(ii)=0;
        RA(ii)=0;
    else
        
        err(ii)=norm(TT-Tref(1:round(dr/drref):end))/...
            norm(Tref(1:round(dr/drref):end));
        RA(ii)=KK/(1-KK*eps)*2*eps;
    end
    
end

figure(1)
grid on
box on
title ('Distribuzione di temperatura')
xlabel('Raggio (cm)')
ylabel('Temperatura (°C)')
set(gca,'fontsize',18)

figure(2)
loglog(drv,err,'-o','linewidth',2)
hold on
loglog(drv,RA,'r-o','linewidth',2)
loglog([drv(2) max(drv)],[1e-4 1e-4],'k--','linewidth',2)
loglog(drv(3),err(3),'ko','markerfacecolor','g','markersize',8)
title('Errore relativo')
xlabel('\deltar (m)')
ylabel('Errore [-]')
set(gca,'fontsize',22)
box on
grid on
legend('Errore rispetto alla soluzione numerica','Limite superiore \epsilon_{RO}',...
    'Tolleranza richiesta','\deltar scelto','location','sw')
xlim([drv(2) max(drv)])