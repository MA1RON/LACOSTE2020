clear
close all
clc

set(0,'defaultaxesfontsize',16)
set(0,'defaultlinelinewidth',2)

%% exercise 2 lab 9 - LaCoSTe

%% data
Lx=1; %m
Ly=1; %m

Th=500; %K
Tc=300; %K
hh=10; %W/m2/K

kk=0.72; %W/m/K
rho=1920; %kg/m3
cp=835; %J/kg/K

%% building the grid

nx=51;
ny=51; 

xvect=linspace(0,Lx,nx);
yvect=linspace(0,Ly,ny);

dx=xvect(2)-xvect(1);
dy=yvect(2)-yvect(1);

[xmat,ymat]=meshgrid(xvect,yvect);
ntot=nx*ny;

%% building the matrix and the known term vector


bb=zeros(ntot,1);

Ac=-2*(1+dx^2/dy^2)*ones(ntot,1);
Ae=ones(ntot,1);
Aw=ones(ntot,1);
An=(dx^2/dy^2)*ones(ntot,1);
As=An;

%set BC

%dirichlet north side
Ae(1:nx)=0;
Aw(1:nx)=0;
An(1:nx)=0;
As(1:nx)=0;

Ac(1:nx)=1;

bb(1:nx)=Th;

%dirichlet west side
Ae(nx+1:nx:end)=0;
Aw(nx+1:nx:end)=0;
An(nx+1:nx:end)=0;
As(nx+1:nx:end)=0;

Ac(nx+1:nx:end)=1;

bb(nx+1:nx:end)=Th;

%dirichlet south side
Ae(nx*(ny-1)+1:end)=0;
Aw(nx*(ny-1)+1:end)=0;
An(nx*(ny-1)+1:end)=0;
As(nx*(ny-1)+1:end)=0;

Ac(nx*(ny-1)+1:end)=1;

bb(nx*(ny-1)+1:end)=Th;

%robin east side
Ae(nx:nx:end)=0;
An(nx:nx:end)=0;
As(nx:nx:end)=0;

Ac(nx:nx:end)=1+(hh*dx)/kk;
Aw(nx:nx:end)=-1;

bb(nx:nx:end)=(hh*dx)/kk*Tc;

%assembling of the matrix

band=[[An(nx+1:end); zeros(nx,1)] [Aw(2:end); 0] Ac [0; Ae(1:end-1)] [zeros(nx,1); As(1:nx*(ny-1))]];
indeces=[-nx -1 0 1 nx];
AA=spdiags(band,indeces,ntot,ntot);

%check the sparsity pattern of the matrix

figure (1)
hold on
spy(AA)

%% solve the steady state problem

TT=AA\bb;

%% arrange the solution in matrix form and plot it

TTmatrix=reshape(TT,[nx,ny]);

figure (2)
hold on
surf(xmat,ymat,TTmatrix','facecolor','interp');
colorbar
xlabel('x [m]');
ylabel('y [m]');
zlabel('z [m]');

%% compute the heat per unit length exchanged with air


qlin=trapz(yvect,hh*(TT(nx:nx:end)-Tc));


disp(['qlin_trapz = ',num2str(qlin),' [W/m]']);


%% TRANSIENT PART!
%to start the transient I need to save the steady state temperature field
%that has been previously calcoulated in order to start the transient.

dt=150; %s
Tf=400; %K
Told=TT;
alpha=kk/rho/cp;

ax=(alpha*dt)/dx^2;
ay=(alpha*dt)/dy^2;

%% compute the matrix coefficients

Act=(1+2*ax*(1+dx^2/dy^2))*ones(ntot,1);
Aet=-ax*ones(ntot,1);
Awt=Aet;
Ant=-ay*ones(ntot,1);
Ast=Ant;

%% fix the boundary conditions

%dirichlet BC north side
Aet(1:nx)=0;
Awt(1:nx)=0;
Ant(1:nx)=0;
Ast(1:nx)=0;

Act(1:nx)=1;

%dirichlet BC west side
Aet(nx+1:nx:end)=0;
Awt(nx+1:nx:end)=0;
Ant(nx+1:nx:end)=0;
Ast(nx+1:nx:end)=0;

Act(nx+1:nx:end)=1;

%dirichlet BC south side
Aet(nx*(ny-1)+1:end)=0;
Awt(nx*(ny-1)+1:end)=0;
Ant(nx*(ny-1)+1:end)=0;
Ast(nx*(ny-1)+1:end)=0;

Act(nx*(ny-1)+1:end)=1;

%robin BC east side
Aet(nx:nx:end)=0;
Ant(nx:nx:end)=0;
Ast(nx:nx:end)=0;

Act(nx:nx:end)=1+(hh*dx)/kk;
Awt(nx:nx:end)=-1;

%% asemble the matrix

band=[[Ant(nx+1:end); zeros(nx,1)] [Awt(2:end); 0] Act [0; Aet(1:end-1)] [zeros(nx,1); Ast(1:nx*(ny-1))]];
indeces=[-nx -1 0 1 nx];
AA=spdiags(band,indeces,ntot,ntot);

figure (3)
hold on
spy(AA);

%% solve the transient

toll=1e-5;
rel_diff=inf;
tempo=0;
ii=1;
ql=[qlin]
tempo=[0]

while rel_diff>toll
    
    
    ii=ii+1;
    tempo(ii)=tempo(ii-1)+dt;
    
    bb=Told;
    bb(1:nx)=Tf;
    bb(nx+1:nx:end)=Tf;
    bb(nx*(ny-1)+1:end)=Tf;
    bb(nx:nx:end)=(hh*dx)/kk*Tc;
    
    %solve
    TT=AA\bb;
    
    ql(ii)=trapz(yvect,hh*(TT(nx:nx:end)-Tc));
    
    rel_diff=abs(ql(ii)-ql(ii-1))/abs(ql(ii));
    
    Tmin=min(TT);
    Told=TT;
    
end

%% plotting

TTmatrix=reshape(TT,[nx,ny]);

figure (4)
hold on
surf(xmat,ymat,TTmatrix','facecolor','interp');
colorbar
xlabel('x [m]');
ylabel('y [m]');
zlabel('z [m]');


figure (5)
hold on
grid minor
xlabel('Time [s]');
ylabel('q_{lin} [W/K]');
plot(tempo,ql);


    
    














