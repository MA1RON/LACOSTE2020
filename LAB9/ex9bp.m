
clear all
close all
clc

%% Lab 9 - Exercise 9

%% Data

ql= 50;      %[W/m]
cond  = 2;       %[W/m/K]
alpha = 1.5e6;   %[m^2/s]
Tair  = 30;      %[°C]
hh    = 100;     %[W/m^2/K]
Lx    = 12e-3;   %[m]    
Ly    = 6e-3;    %[m]
Lyh   = 2e-3;    %[m]

%% Grid Building
% Scelgo una porzione della piastra che includa un singolo riscaldatore
% sul contorno sinistro del dominio computazionale

dx = 1e-3;
dy = dx/2;

xvec = [0:dx:Lx]';
yvec = [0:dy:Ly]';
nx = length(xvec);
ny = length(yvec);
nh = round(Lyh/dy)+1;        %index of the heater position along the west side

[xmat,ymat] = meshgrid(xvec,yvec);
ntot = nx*ny;

%% Matrix and Known Vector Building

%Sparse Matrix and Known Vector Allocation
AA = sparse([],[],[],ntot,ntot,5*ntot);
rhs = zeros(ntot,1);

%Internal points
for icenter = 2:nx-1
    for jcenter = 2:ny-1
        kcenter = (jcenter-1)*nx+icenter;
        knorth  = kcenter+nx;
        ksouth  = kcenter-nx;
        keast   = kcenter+1;
        kwest   = kcenter-1;
        
        AA(kcenter,kcenter) = -2*cond*(1/dx^2+1/dy^2);
        AA(kcenter,keast)   = cond/(dx^2);
        AA(kcenter,kwest)   = cond/(dx^2);
        AA(kcenter,knorth)  = cond/(dy^2);
        AA(kcenter,ksouth)  = cond/(dy^2);
    end
end

%Boundary Conditions

%South - Adiabatic
jcenter = 1;
for icenter = 2:nx-1
    kcenter = (jcenter-1)*nx+icenter;
    knorth  = kcenter+nx;
    keast   = kcenter+1;
    kwest   = kcenter-1;

    AA(kcenter,kcenter) = cond*(-1/(dy^2)-1/(dx^2));
    AA(kcenter,keast)   = cond/(2*dx^2);
    AA(kcenter,kwest)   = cond/(2*dx^2);
    AA(kcenter,knorth)  = cond/(dy^2);
    
    rhs(kcenter) = 0;
end

%East - Adiabatic
icenter = nx;
for jcenter = 2:ny-1
    kcenter = (jcenter-1)*nx+icenter;
    knorth  = kcenter+nx;
    ksouth  = kcenter-nx;
    kwest   = kcenter-1;
    
    AA(kcenter,kcenter)  = cond*(-1/(dy^2)-1/(dx^2));
    AA(kcenter,ksouth)   = cond/(2*dy^2);
    AA(kcenter,kwest)    = cond/(dx^2);
    AA(kcenter,knorth)   = cond/(2*dy^2);
    
    rhs(kcenter) = 0;
end

%West - Adiabatic + Heat Supply
icenter = 1;
for jcenter = 2:ny-1
    kcenter = (jcenter-1)*nx+icenter;
    knorth  = kcenter+nx;
    ksouth  = kcenter-nx;
    keast   = kcenter+1;
    
    AA(kcenter,kcenter)  = cond*(-1/(dy^2)-1/(dx^2));
    AA(kcenter,ksouth)   = cond/(2*dy^2);
    AA(kcenter,keast)    = cond/(dx^2);
    AA(kcenter,knorth)   = cond/(2*dy^2);
    
    rhs(kcenter) = 0;
    
    if jcenter == nh
        rhs(kcenter) = -ql/(2*dx*dy);
    end
end

%North - Robin
jcenter = ny;
for icenter = 2:nx-1
    kcenter = (jcenter-1)*nx+icenter;
    ksouth  = kcenter-nx;
    keast   = kcenter+1;
    kwest   = kcenter-1;

    AA(kcenter,kcenter) = cond*(-1/(dy^2)-1/(dx^2))-hh/dy;
    AA(kcenter,keast)   = cond/(2*dx^2);
    AA(kcenter,kwest)   = cond/(2*dx^2);
    AA(kcenter,ksouth)  = cond/(dy^2);
    
    rhs(kcenter) = -hh*Tair/dy;
end

%South-West Corner
icenter = 1;
jcenter = 1;
kcenter = (jcenter-1)*nx+icenter;
keast   = kcenter+1;
knorth  = kcenter+nx;

AA(kcenter,kcenter) = cond*(-1/(dy^2)-1/(dx^2));
AA(kcenter,keast)   = cond/(dx^2);
AA(kcenter,knorth)  = cond/(dy^2);

rhs(kcenter) = 0;

%South-East Corner
icenter = nx;
jcenter = 1;
kcenter = (jcenter-1)*nx+icenter;
kwest   = kcenter-1;
knorth  = kcenter+nx;

AA(kcenter,kcenter) = cond*(-1/(dy^2)-1/(dx^2));
AA(kcenter,kwest)   = cond/(dx^2);
AA(kcenter,knorth)  = cond/(dy^2);

rhs(kcenter) = 0;

%North-East Corner
icenter = nx;
jcenter = ny;
kcenter = (jcenter-1)*nx+icenter;
kwest   = kcenter-1;
ksouth  = kcenter-nx;

AA(kcenter,kcenter) = cond*(-1/(dy^2)-1/(dx^2))-hh/dy;
AA(kcenter,kwest)   = cond/(dx^2);
AA(kcenter,ksouth)  = cond/(dy^2);

rhs(kcenter) = -hh*Tair/dy;

%North-West Corner
icenter = 1;
jcenter = ny;
kcenter = (jcenter-1)*nx+icenter;
keast   = kcenter+1;
ksouth  = kcenter-nx;

AA(kcenter,kcenter) = cond*(-1/(dy^2)-1/(dx^2))-hh/dy;
AA(kcenter,keast)   = cond/(dx^2);
AA(kcenter,ksouth)  = cond/(dy^2);

rhs(kcenter) = -hh*Tair/dy;


%% Steady state Solution and Plot

DeltaTT0 = AA\rhs - Tair;
TTmatrix = reshape(DeltaTT0,nx,ny);

figure (1)
hold on
surf(xmat*1e3,ymat*1e3,TTmatrix')
colorbar
xlabel('x [mm]');
ylabel('y [mm]');
title('T-Tair [°C]');

figure(2)
hold on
contourf(xmat*1e3,ymat*1e3,TTmatrix')
colorbar
xlabel('x [mm]');
ylabel('y [mm]');
title('T-Tair [°C]');

%% Transient State

deltat = 2; %[s]
BB = eye(length(AA))-AA*alpha*deltat/cond;  %nella definizione di AA è già contenuta la conducibilità!!
qv = zeros(length(BB),1);

%South - Adiabatic
jcenter = 1;
for icenter = 2:nx-1
    kcenter = (jcenter-1)*nx+icenter;
    knorth  = kcenter+nx;
    keast   = kcenter+1;
    kwest   = kcenter-1;

    BB(kcenter,kcenter) = cond*(-1/(dy^2)-1/(dx^2));
    BB(kcenter,keast)   = cond/(2*dx^2);
    BB(kcenter,kwest)   = cond/(2*dx^2);
    BB(kcenter,knorth)  = cond/(dy^2);
end

%East - Adiabatic
icenter = nx;
for jcenter = 2:ny-1
    kcenter = (jcenter-1)*nx+icenter;
    knorth  = kcenter+nx;
    ksouth  = kcenter-nx;
    kwest   = kcenter-1;
    
    BB(kcenter,kcenter)  = cond*(-1/(dy^2)-1/(dx^2));
    BB(kcenter,ksouth)   = cond/(2*dy^2);
    BB(kcenter,kwest)    = cond/(dx^2);
    BB(kcenter,knorth)   = cond/(2*dy^2);
end

%West - Adiabatic + Heat Supply
icenter = 1;
for jcenter = 2:ny-1
    kcenter = (jcenter-1)*nx+icenter;
    knorth  = kcenter+nx;
    ksouth  = kcenter-nx;
    keast   = kcenter+1;
    
    BB(kcenter,kcenter)  = cond*(-1/(dy^2)-1/(dx^2));
    BB(kcenter,ksouth)   = cond/(2*dy^2);
    BB(kcenter,keast)    = cond/(dx^2);
    BB(kcenter,knorth)   = cond/(2*dy^2);
    if jcenter == nh
            qv(kcenter) = ql/(2*dx*dy);
    end
end

%North - Robin
jcenter = ny;
for icenter = 2:nx-1
    kcenter = (jcenter-1)*nx+icenter;
    ksouth  = kcenter-nx;
    keast   = kcenter+1;
    kwest   = kcenter-1;

    BB(kcenter,kcenter) = cond*(-1/(dy^2)-1/(dx^2))-hh/dy;
    BB(kcenter,keast)   = cond/(2*dx^2);
    BB(kcenter,kwest)   = cond/(2*dx^2);
    BB(kcenter,ksouth)  = cond/(dy^2);
end

%South-West Corner
icenter = 1;
jcenter = 1;
kcenter = (jcenter-1)*nx+icenter;
keast   = kcenter+1;
knorth  = kcenter+nx;

BB(kcenter,kcenter) = cond*(-1/(dy^2)-1/(dx^2));
BB(kcenter,keast)   = cond/(dx^2);
BB(kcenter,knorth)  = cond/(dy^2);

%South-East Corner
icenter = nx;
jcenter = 1;
kcenter = (jcenter-1)*nx+icenter;
kwest   = kcenter-1;
knorth  = kcenter+nx;

BB(kcenter,kcenter) = cond*(-1/(dy^2)-1/(dx^2));
BB(kcenter,kwest)   = cond/(dx^2);
BB(kcenter,knorth)  = cond/(dy^2);

%North-East Corner
icenter = nx;
jcenter = ny;
kcenter = (jcenter-1)*nx+icenter;
kwest   = kcenter-1;
ksouth  = kcenter-nx;

BB(kcenter,kcenter) = cond*(-1/(dy^2)-1/(dx^2))-hh/dy;
BB(kcenter,kwest)   = cond/(dx^2);
BB(kcenter,ksouth)  = cond/(dy^2);

%North-West Corner
icenter = 1;
jcenter = ny;
kcenter = (jcenter-1)*nx+icenter;
keast   = kcenter+1;
ksouth  = kcenter-nx;

BB(kcenter,kcenter) = cond*(-1/(dy^2)-1/(dx^2))-hh/dy;
BB(kcenter,keast)   = cond/(dx^2);
BB(kcenter,ksouth)  = cond/(dy^2);


DeltaTT0_av = trapz(xvec,DeltaTT0(end-nx+1:end))/Tair %Average Temperature Difference on the North surface
DeltaTT_av = 0.01*DeltaTT0_av;
time=0;
TTold =  Tair*ones(size(DeltaTT0));
counter=0;
while DeltaTT_av < 0.95*DeltaTT0_av
    counter=counter+1;
    time = time + deltat
    rhs  = TTold + qv*deltat*alpha/cond;
    
    %South - Adiabatic
        rhs(2:nx-1) = 0;

    %East - Adiabatic
    icenter = nx;
    for jcenter = 2:ny-1
        kcenter = (jcenter-1)*nx+icenter;
        rhs(kcenter) = 0;
    end

    %West - Adiabatic + Heat Supply
    icenter = 1;
    for jcenter = 2:ny-1
        kcenter = (jcenter-1)*nx+icenter;
        if jcenter ~= nh
            rhs(kcenter) = 0;
        end
    end

    %North - Robin
    jcenter = ny;
    for icenter = 2:nx-1
        kcenter = (jcenter-1)*nx+icenter;
        rhs(kcenter) = -hh*Tair/dy;
    end

    %South-West Corner
    rhs(1) = 0;

    %South-East Corner
    rhs(nx) = 0;

    %North-East Corner
    rhs(ny*nx) = -hh*Tair/dy;

    %North-West Corner
    rhs((ny-1)*nx+1) = -hh*Tair/dy;
    
    TTnew = BB\rhs;
    Temp(counter)=TTnew(end-nx+1);
    tempo(counter)=time;
    DeltaTT_av = abs(trapz(xvec,TTnew(end-nx+1:end)-Tair))/Tair
    TTmatrix = reshape(TTnew,nx,ny);

    figure(1)
    surf(xmat*1e3,ymat*1e3,TTmatrix')
    pause(0.2)
    TTold = TTnew;
end

request_time = time;









