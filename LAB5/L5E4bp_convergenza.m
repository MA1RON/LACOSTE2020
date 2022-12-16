clc
clear 
close all

%% Assegno i dati
LL=50e-2; %m
bb=10e-2; %m
th=2e-2; %m
kk=350; %W/(m*K)
hh=500; %W/(m^2*K)
The=4.5; %K
II=7e3; %A
rho_el0=1.75e-8; %ohm*m
L0=56e-2; %m

rhoel=@(x) rho_el0*(cos(pi*x./L0)).^2; %ohm*m
As=2*(bb+th)*LL; %m^2
SS=th*bb; %m^2
VV=SS*LL; %m^3

%% Studio di convergenza
dxv=[1e-6 5e-6 1e-5 5e-5 1e-4 5e-4 1e-3];
err=size(dxv);
RA=err;

for ii=1:length(dxv)

    dx=dxv(ii);
    xx=(0:dx:LL/2)';
    nn=length(xx);
    
    qvol=II^2*LL*rhoel(xx)./(SS*VV);
    CC=hh*As/VV/kk;
    
    sub_diag=ones(nn,1);
    main=-(2+CC*dx^2)*ones(nn,1);
    super_diag=sub_diag;
    
    BB=[sub_diag,main,super_diag];
    
    AA=spdiags(BB,-1:1,nn,nn);
    
    bb=-qvol*dx^2./kk-CC*dx^2*The;
    
    AA(1,1)=1;
    AA(1,2)=-1;
    bb(1)=0;
    
    AA(end,end-1)=0;
    AA(end,end)=1;
    bb(end)=The;
    
    TT=AA\bb;
    
    KK=condest(AA);
    
    figure(1)
    plot(xx,TT,'displayname',['dx=' num2str(dxv(ii)) ' m'],...
        'linewidth',2)
    hold on
    legend('-dynamiclegend','location','sw')
    
    if ii==1
        
        dxref=dx;
        Tref=TT;
        
        err(ii)=0;
        RA(ii)=0;
    else
        
        err(ii)=norm(TT-Tref(1:round(dx/dxref):end))/...
            norm(Tref(1:round(dx/dxref):end));
        RA(ii)=KK/(1-KK*eps)*2*eps;
    end
    

end

figure(1)
title('Distribuzione di temperatura')
xlabel('Lunghezza (m)')
ylabel('Temperatura (K)')
set(gca,'fontsize',18)
box on
grid on

figure(2)
loglog(dxv,err,'-o','linewidth',2)
hold on
loglog(dxv,RA,'r-o','linewidth',2)
loglog([dxv(2) max(dxv)],[1e-4 1e-4],'k--','linewidth',2)
loglog(dxv(5),err(5),'ko','markerfacecolor','g','markersize',8)
title('Errore relativo flusso termico')
xlabel('\deltax (m)')
ylabel('Errore [-]')
set(gca,'fontsize',18)
box on
grid on
legend('Errore rispetto alla soluzione numerica','Limite superiore \epsilon_{RO}',...
    'Tolleranza richiesta','\deltax scelto','location','sw')