clear 
close all
clc

set(0,'defaultaxesfontsize',30)
set(0,'defaultlinelinewidth',2)

%% lacoste9_3
%exercise 3 lab 9 lacoste

%% data
kk=1.5;
Lx=0.4;
Ly=0.6;
Tf=200; %°C
hh=50;
Ta=30;

%% building the grid

nx=101;
ny=101; 

xvect=linspace(0,Lx,nx);
yvect=linspace(0,Ly,ny);

dx=xvect(2)-xvect(1);
dy=yvect(2)-yvect(1);

[xmat,ymat]=meshgrid(xvect,yvect);
ntot=nx*ny;


%% building the matrix diagonals

bb=zeros(ntot,1);

Ac=-2*(1+dx^2/dy^2)*ones(ntot,1);
Ae=ones(ntot,1);
Aw=Ae;
An=(dx^2/dy^2)*ones(ntot,1);
As=An;

%% now fix the BC

%dirichlet north side
Ae(1:nx)=0;
Aw(1:nx)=0;
An(1:nx)=0;
As(1:nx)=0;

Ac(1:nx)=1;

bb(1:nx)=Tf;

%homogeneus neumann west side
Ae(nx+1:nx:end)=1;
Ac(nx+1:nx:end)=-1;
Aw(nx+1:nx:end)=0;
An(nx+1:nx:end)=0;
As(nx+1:nx:end)=0;

bb(nx+1:nx:end)=0;

%dirichlet south side
Ae(nx*(ny-1)+1:end)=0;
Aw(nx*(ny-1)+1:end)=0;
An(nx*(ny-1)+1:end)=0;
As(nx*(ny-1)+1:end)=0;

Ac(nx*(ny-1)+1:end)=1;

bb(nx*(ny-1)+1:end)=Tf;

%robin east side
Ae(nx:nx:nx*(ny-1))=0;
An(nx:nx:nx*(ny-1))=0;
As(nx:nx:nx*(ny-1))=0;

Ac(nx:nx:nx*(ny-1))=1+(hh*dx)/kk;
Aw(nx:nx:nx*(ny-1))=-1;

bb(nx:nx:nx*(ny-1))=(hh*dx)/kk*Ta;

%assembling of the matrix

band=[[An(nx+1:end); zeros(nx,1)] [Aw(2:end); 0] Ac [0; Ae(1:end-1)] [zeros(nx,1); As(1:nx*(ny-1))]];
indeces=[-nx -1 0 1 nx];
AA=spdiags(band,indeces,ntot,ntot);


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

figure (3)
hold on
contourf(xmat,ymat,TTmatrix',30);
colorbar
xlabel('x [m]');
ylabel('y [m]');

%% compute the heat per unit length exchanged with air


qlin=trapz(yvect,hh*(TT(nx:nx:end)-Ta));


disp(['qlin_trapz = ',num2str(qlin),' [W/m]']);


