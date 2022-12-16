%% LaCoSTe 2020/2021 
%% Lab 9 Ex 6 - Soluzioni

clear all
close all
clc

cond=50;
w=6e-3;
L=48e-3;

h=[10,100,500,1000];
Ta=30;
Tb=100;

dx=5e-4;
dy=dx;

xvec=(0:dx:L)';
yvec=(0:dy:w/2)';
[xmat,ymat]=meshgrid(xvec,yvec);

nx=length(xvec);
ny=length(yvec);

ntot=nx*ny;
A=sparse([],[],[],ntot,ntot);
rhs=zeros(ntot,1);

%% Matrice dei coefficienti
for icenter=2:nx-1
   for jcenter=2:ny-1
    kcenter=(jcenter-1)*nx+icenter;
    knorth=kcenter+nx;
    ksouth=kcenter-nx;
    keast=kcenter+1;
    kwest=kcenter-1;

    A(kcenter,kcenter)=2*cond*(1/dx^2+1/dy^2);
    A(kcenter,knorth)=-cond/dy^2;
    A(kcenter,ksouth)=-cond/dy^2;
    A(kcenter,keast)=-cond/dx^2;
    A(kcenter,kwest)=-cond/dx^2;
   end
end

%% Condizioni al contorno
% Dirichlet
icenter=1;
for jcenter=1:ny
    kcenter=(jcenter-1)*nx+icenter;
    A(kcenter,kcenter)=1.0;
    rhs(kcenter)=Tb;
end
% Neumann adiabatico
jcenter=1;
for icenter=2:nx-1
    kcenter=(jcenter-1)*nx+icenter;
    knorth=kcenter+nx;
    kwest=kcenter-1;
    keast=kcenter+1;
    A(kcenter,kcenter)=2*cond/dy^2+2*cond/dx^2;
    A(kcenter,knorth)=-2*cond/dy^2;
    A(kcenter,kwest)=-cond/dx^2;
    A(kcenter,keast)=-cond/dx^2;
    rhs(kcenter)=0;
end
% Neumann adiabatico
icenter=nx;
for jcenter=2:ny-1
    kcenter=(jcenter-1)*nx+icenter;
    knorth=kcenter+nx;
    kwest=kcenter-1;
    ksouth=kcenter-nx;
    A(kcenter,kcenter)=2*cond/dy^2+2*cond/dx^2;
    A(kcenter,knorth)=-cond/dy^2;
    A(kcenter,kwest)=-2*cond/dx^2;
    A(kcenter,ksouth)=-cond/dy^2;
    rhs(kcenter)=0;
end

kcenter=nx;
knorth=kcenter+nx;
kwest=kcenter-1;
A(kcenter,kcenter)=cond/dy^2+cond/dx^2;
A(kcenter,knorth)=-cond/dy^2;
A(kcenter,kwest)=-cond/dx^2;
rhs(kcenter)=0;

for i=1:length(h)
    % Robin
    jcenter=ny;
    for icenter=2:nx-1
        kcenter=(jcenter-1)*nx+icenter;
        ksouth=kcenter-nx;
        kwest=kcenter-1;
        keast=kcenter+1;
        A(kcenter,kcenter)=2*(cond/dy^2+cond/dx^2+h(i)/dy);
        A(kcenter,keast)=-cond/dx^2;
        A(kcenter,kwest)=-cond/dx^2;
        A(kcenter,ksouth)=-2*cond/dy^2;
        rhs(kcenter)=2*h(i)*Ta/dy;
    end

    kcenter=nx*ny;
    ksouth=kcenter-nx;
    kwest=kcenter-1;
    A(kcenter,kcenter)=-(cond/dy^2+cond/dx^2+h(i)/dx);
    A(kcenter,ksouth)=cond/dy^2;
    A(kcenter,kwest)=cond/dx^2;
    rhs(kcenter)=-h(i)/dx*Ta;
    
    if h(i)==500 %Controllo le condizioni al contorno (per h=500W/m2K)
        figure(1)
        rrr=reshape(rhs,nx,ny);
        surf(xmat,ymat,rrr')
        grid minor
        title('Controllo condizioni al contorno')
    end
    
    T=A\rhs;
    figure(3)
    grid minor
    plot(xvec,T(1:nx),'Linewidth',2)
    hold on
    xlabel('x [m]')
    ylabel('T [°C]')

    if h(i)==500 %Plotto la distribuzione di temperatura (per h=500W/m2K)
        figure(2)
        T_reshape=reshape(T,nx,ny);
        surf(xmat,ymat,T_reshape')
        hold on
        grid minor
        title('Distribuzione di temperatura')
        xlabel('x [m]')
        ylabel('y [m]')
        zlabel('T [°C]')
        
        %Calcolo il calore dissipato per convezione, considerando anche la
        %seconda faccia dell'aletta.
        q=2*trapz(xvec,h(i)*(T_reshape(:,end)-Ta)); %[W/m]
    end
end

figure(3)
title('Profilo di temperatura lungo asse di simmetria aletta')
legend('h=10 W/m^2K','h=100 W/m^2K','h=500 W/m^2K','h=1000 W/m^2K','Location','Southwest')


