clear all
close all
clc


%% Data

Lx  = 30e-3;   %[m]
Ly1 = 50e-3;   %[m]
Ly2 = 20e-3;   %[m]
T1  = 100;     %[°C]
T2  = 25;      %[°C]
cond  = 20;    %[W/m/K]


%% Grid Building

nx   = 31;
ny1  = 51;
ny2  = 21;

xvec = linspace(0,Lx,nx);
yvec = linspace(0,Ly1,ny1);
dx = xvec(2)-xvec(1);
dy = yvec(2)-yvec(1);

[xmat,ymat] = meshgrid(xvec,yvec);
ntot = nx*(ny2-1)+(ny1-ny2)*nx/2+nx;

% temp1(1:nx,1:ny1) = 0;
% 
% for ii = 2:nx
%     for jj = ny1-ii+2:ny1
%             temp1(ii,jj) = NaN;
%     end
% end
% 
% surf(xmat,ymat,temp1');

%% Matrix and Known Vector Building

%Sparse Matrix and Known Vector Allocation
AA = sparse([],[],[],ntot,ntot,5*ntot);
rhs = zeros(ntot,1);

%Internal Points - Rectangle
for icenter = 2:nx-1
    for jcenter = 2:ny2
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

%Internal Points - Triangle
nbase = nx*ny2;               %somma degli elementi che precedono ciascuna riga
ni = nx-1;                    %numero degli elementi sulla riga in esame
for jcenter = ny2+1:ny1-2  
    for icenter = 2:ni-1      %sottraggo il bordo del triangolo
        kcenter = nbase+icenter;
        knorth  = kcenter+ni;
        ksouth  = kcenter-(ni+1);
        keast   = kcenter+1;
        kwest   = kcenter-1;

        AA(kcenter,kcenter) = -2*cond*(1/dx^2+1/dy^2);
        AA(kcenter,keast)   = cond/(dx^2);
        AA(kcenter,kwest)   = cond/(dx^2);
        AA(kcenter,knorth)  = cond/(dy^2);
        AA(kcenter,ksouth)  = cond/(dy^2);  
    end
    nbase = nbase+ni;
    ni = ni-1;
end

%nbase_uppercorner = nbase;
    

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

%East - Dirichelet
icenter = nx;
for jcenter = 1:ny2
    kcenter = (jcenter-1)*nx+icenter;

    AA(kcenter,kcenter) = 1;
    rhs(kcenter) = T2;
end

%West - Dirichelet
icenter = 1;
for jcenter = 1:ny2
    kcenter = (jcenter-1)*nx+icenter;

    AA(kcenter,kcenter)  = 1;
    rhs(kcenter) = T1;
end

%Diagonal - Adiabatic
nbase = nx*ny2;
ni = nx-1;
for jcenter = ny2+1:ny1-1
    
    kcenter = nbase+ni;
    ksouth  = kcenter-(ni+1);
    kwest   = kcenter-1;

    AA(kcenter,kcenter) = -cond*(1/dx^2+1/dy^2);
    AA(kcenter,kwest)   = cond/(dx^2);
    AA(kcenter,ksouth)  = cond/(dy^2);  
    
    rhs(kcenter) = 0;
    
    nbase = nbase+ni;
    ni = ni-1;

end


%West Triangle - Dirichelet
nbase = nx*ny2;
ni = nx-1;
for jcenter = ny2+1:ny1
    icenter = 1;
    kcenter = nbase+icenter;

    AA(kcenter,kcenter)  = 1;
    rhs(kcenter) = T1;
    
    nbase = nbase+ni;
    ni = ni - 1;
end


% %Triangular Upper Corner 
% icenter = 2;
% jcenter = ny2-2;
% ni=3;
% 
% kcenter = nbase_uppercorner+icenter;
% knorth  = kcenter+ni;
% ksouth  = kcenter-(ni+1);
% keast   = kcenter+1;
% kwest   = kcenter-1;
% 
% AA(kcenter,kcenter) = -2*cond*(1/dx^2+1/dy^2);
% AA(kcenter,keast)   = cond/(dx^2);
% AA(kcenter,knorth)  = cond/(dy^2);
% AA(kcenter,ksouth)  = cond/(dy^2); 
% 
% rhs(kcenter) = -T1/(dx^2);



%% Temperature Distribution and Plot

TT = AA\rhs;

nbase = nx*ny2;
TTmatrix1 = reshape(TT(1:nbase),nx,ny2);


ni = nx-1;                    
for j = 1:ny1-ny2
    TTmatrix2(j,:) = [TT(nbase+1:nbase+ni)',NaN*ones(1,nx-ni)];
    
    nbase = nbase+ni;
    ni = ni-1;
end
TTmatrix = [TTmatrix1,TTmatrix2'];

figure (1)
hold on
surf(xmat,ymat,TTmatrix')
colorbar
xlabel('x [m]');
ylabel('y [m]');

figure(2)
hold on
contourf(xmat,ymat,TTmatrix')
colorbar
xlabel('x [m]');
ylabel('y [m]');

%% Heat Transfer Evaluation

%East Side
icenter = nx;
qq_out = 0;
for jcenter = 1:ny2
    kcenter = (jcenter-1)*nx+icenter;
    kwest   = kcenter-1;
    
    Tw = TT(kwest);
    Tc = TT(kcenter);
    
    if jcenter == 1 || jcenter == ny2
        qq_out = qq_out + cond*(Tw-Tc)*dy/(2*dx);
    else
        qq_out = qq_out + cond*(Tw-Tc)*dy/dx;
    end
end


%West Side
icenter = 1;
qq_in = 0;
for jcenter = 1:ny2
    kcenter = (jcenter-1)*nx+icenter;
    keast   = kcenter+1;
    
    Te = TT(keast);
    Tc = TT(kcenter);
    
    if jcenter == 1 || jcenter == ny2
        qq_in = qq_in + cond*(Tc-Te)*dy/(2*dx);
    else
        qq_in = qq_in + cond*(Tc-Te)*dy/dx;
    end
end

nbase = nx*ny2;
ni = nx-1;
for jcenter = ny2+1:ny1-1
    icenter = 1;
    kcenter = nbase+icenter;
    keast   = kcenter+1;

    Te = TT(keast);
    Tc = TT(kcenter);
    
    qq_in = qq_in + cond*(Tc-Te)*dy/dx;

    nbase = nbase+ni;
    ni = ni - 1;
end






