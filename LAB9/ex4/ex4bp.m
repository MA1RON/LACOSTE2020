clear 
close all
clc

set(0,'defaultaxesfontsize',14)
set(0,'defaultlinelinewidth',2)

%% lacoste9_4
%exercise 4 lab 9 lacoste

%% data
LL=27e-3;
HH=12e-3;
ks=5;
kc=50;
Ta=20; %°C
hh=500;
qv0=1e7;
Tmax=85; %°C

%% generation of the grid
%I generate the grid of the overall domain 
%--> COULD THE COMPUTATIONAL DOMAIN BE REDUCED?

nx=51;
ny=51; 

xvect=linspace(0,LL,nx);
yvect=linspace(0,HH,ny);

dx=xvect(2)-xvect(1);
dy=yvect(2)-yvect(1);

[xmat,ymat]=meshgrid(xvect,yvect);
ntot=nx*ny;

kfunc=@(x,y) ks+(kc-ks).*(x>=9e-3 & x<=18e-3 & y>=9e-3);
qvfunc=@(x,y) 0+qv0.*(x>=9e-3 & x<=18e-3 & y>=9e-3);

kmat=kfunc(xmat,ymat);
qvmat=qvfunc(xmat,ymat);

figure (1)
hold on
title('Thermal conductivity');
xlabel('x [m]');
ylabel('y [m]');
zlabel('k [W/m/K]');
surf(xmat,ymat,kmat);
colorbar

figure (2)
hold on
title('Volumetric heat generation');
xlabel('x [m]');
ylabel('y [m]');
zlabel('qv [W/m^3]');
surf(xmat,ymat,qvmat);
colorbar

%voglio kvect --> kk
kk=flip(reshape(kmat',[ntot,1]));
qv=flip(reshape(qvmat',[ntot,1]));

%% find a way to solve the equation

%known term vector
bb=-(qv*dx^2)./kk;

Ac=-2*(1+dx^2/dy^2)*ones(ntot,1);
Ae=ones(ntot,1);
Aw=Ae;
An=(dx^2/dy^2)*ones(ntot,1);
As=An;

% set boundary conditions

%neumann west side
Ae(1:nx:end)=1;
Ac(1:nx:end)=-1;
bb(1:nx:end)=0;

Aw(1:nx:end)=0;
As(1:nx:end)=0;
An(1:nx:end)=0;

%neumann south side
An(nx*(ny-1)+2:end)=1;
Ac(nx*(ny-1)+2:end)=-1;
bb(nx*(ny-1)+2:end)=0;

Aw(nx*(ny-1)+2:end)=0;
As(nx*(ny-1)+2:end)=0;
Ae(nx*(ny-1)+2:end)=0;

%neumann east side
Ac(nx:nx:(ny-1)*nx)=1;
Aw(nx:nx:(ny-1)*nx)=-1;
bb(nx:nx:(ny-1)*nx)=0;

An(nx:nx:(ny-1)*nx)=0;
As(nx:nx:(ny-1)*nx)=0;
Ae(nx:nx:(ny-1)*nx)=0;

%robin north side
Ac(2:nx-1)=1+(hh*dy)./kk(2:nx-1)
As(2:nx-1)=-1;
bb(2:nx-1)=(hh*dy)./kk(2:nx-1)*Ta;

An(2:nx-1)=0;
Ae(2:nx-1)=0;
Aw(2:nx-1)=0;

%now I need to impose the interface conditions

%interface 1 (west)
Ac(ceil(nx/3)+1+nx:nx:ceil(nx/3)+1+(ceil(ny/4)-1)*nx)=1+kc/ks;
Aw(ceil(nx/3)+1+nx:nx:ceil(nx/3)+1+(ceil(ny/4)-1)*nx)=-1;
Ae(ceil(nx/3)+1+nx:nx:ceil(nx/3)+1+(ceil(ny/4)-1)*nx)=-kc/ks
bb(ceil(nx/3)+1+nx:nx:ceil(nx/3)+1+(ceil(ny/4)-1)*nx)=0;

As(ceil(nx/3)+1+nx:nx:ceil(nx/3)+1+(ceil(ny/4)-1)*nx)=0;
An(ceil(nx/3)+1+nx:nx:ceil(nx/3)+1+(ceil(ny/4)-1)*nx)=0;

%interface 2 (south)

puntoA=ceil(nx/3)+1+(ceil(ny/4)-1)*nx+1;
puntoB=puntoA+ceil(nx/3)-1;
Ac(puntoA:puntoB)=1+kc/ks;
As(puntoA:puntoB)=-1;
An(puntoA:puntoB)=-kc/ks;
bb(puntoA:puntoB)=0;

Ae(puntoA:puntoB)=0;
Aw(puntoA:puntoB)=0;

%interface 3 (east)
puntoA=2*ceil(nx/3)-1+nx;
PuntoB=puntoA+(ceil(ny/4)-3)*nx;
Ac(puntoA:nx:puntoB)=1+ks/kc;
Aw(puntoA:nx:puntoB)=-1;
Ae(puntoA:nx:puntoB)=-ks/kc;
bb(puntoA:nx:puntoB)=0;

An(puntoA:nx:puntoB)=0;
As(puntoA:nx:puntoB)=0;

%% matrix assembling
band=[[An(nx+1:end); zeros(nx,1)] [Aw(2:end); 0] Ac [0; Ae(1:end-1)] [zeros(nx,1); As(1:nx*(ny-1))]];
indeces=[-nx -1 0 1 nx];
AA=spdiags(band,indeces,ntot,ntot);

%% solve the steady state problem

TT=AA\bb;

%% arrange the solution in matrix form and plot it

TTmatrix=reshape(flip(TT),[nx,ny]);

figure (3)
hold on
surf(xmat,ymat,TTmatrix','facecolor','interp');
colorbar
xlabel('x [m]');
ylabel('y [m]');
zlabel('z [m]');

figure (4)
hold on
contourf(xmat,ymat,TTmatrix',30);
colorbar
xlabel('x [m]');
ylabel('y [m]');











