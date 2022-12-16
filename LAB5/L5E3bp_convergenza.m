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

Tbottom=70; %K
Tair=20; %K
hh=15; %W/(m^2*K)

%Resistenze termiche
As=0.6*0.6; %m^2
Rss=tss/kss/As; %K/W
Rins=tins/kins/As; %K/W
R=Rss+Rins; %K/W
U=1/R/As; %W/(m^2*K)

%% Studio di convergenza
dxv=[1.5e-7 1.5e-6 1.5e-5 1.5e-4];
err=size(dxv);
err_q=size(dxv);
cond=size(dxv);

for ii=1:length(dxv)
    
    dx=dxv(ii);
    
    xss=(0:dx:tss)';
    NN=length(xss);
    
    xins=(tss+dx:dx:tins+tss)';
    
    xx=[xss;xins];
    nn=length(xx);
    
    kk=kss*(xx<=tss)+kins*(xx>tss);
    
    sub_diag=ones(nn,1);
    main=-2*ones(nn,1);
    super_diag=sub_diag;
    
    BB=[sub_diag,main,super_diag];
    
    AA=spdiags(BB,-1:1,nn,nn);
    
    bb=zeros(nn,1);
    
    AA(1,1)=1;
    AA(1,2)=0;
    bb(1)=Tbottom;
    
    AA(NN,NN-1)=kss/kins;
    AA(NN,NN)=-1-kss/kins;
    AA(NN,NN+1)=1;
    bb(NN)=0;
    
    AA(end,end-1)=-1;
    AA(end,end)=1+hh*dx/kins;
    bb(end)=hh*dx*Tair/kins;
    
    TT=AA\bb;
    
    q1=abs(hh*(TT(end)-Tair));
    q2=abs(U*(Tbottom-TT(end)));
    
    err_q(ii)=abs(q1-q2);
    
    KK=condest(AA);
    
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
loglog(dxv,err,'-o','linewidth',2)
hold on
loglog(dxv,RA,'r-o','linewidth',2)
loglog([dxv(2) dxv(end)],[1e-4 1e-4],'k--','linewidth',2)
loglog(dxv(3),err(3),'ko','markerfacecolor','g','markersize',8)
title('Errore rispetto alla soluzione più raffinata')
xlabel('\deltax (m)')
ylabel('Error (-)')
set(gca,'fontsize',16)
box on
grid on
legend('Errore rispetto alla soluzione numerica','Limite superiore \epsilon_{RO}',...
    'Tolleranza richiesta','\deltax scelto','location','sw')
xlim([dxv(2) dxv(end)])

figure(2)
loglog(dxv,err_q,'r--o','linewidth',2)
title('Errore relativo flusso termico')
xlabel('\deltax (m)')
ylabel('Errore (-)')
set(gca,'fontsize',18)
grid on