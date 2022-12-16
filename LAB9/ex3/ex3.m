clc; clear; close all;
set(0,'defaultaxesfontsize',15)
set(0,'defaultlinelinewidth',2)
%% dati
Lx = .4;
Ly = .6;
cond = 1.5;
Tnorth = 200;
Tsouth = Tnorth;
qqq = 0;
qqwest = 0;
hheast = 50;
Teast = 30;

%% griglia
nx = 101;
ny = 101;
ntot = nx*ny;
xvet = linspace(0,Lx,nx);
yvet = linspace(0,Ly,ny);
dx = xvet(2)-xvet(1);
dy = yvet(2)-yvet(1);

[xmat, ymat] = meshgrid(xvet,yvet);

%% creo il sistema
AA = sparse([],[],[],ntot,ntot,5*ntot);
bb = zeros(ntot,1);

% --- center ---
for jx = 2:nx-1
    for jy = 2:ny-1
        kk = jx+(jy-1)*nx;
        AA(kk,kk) = -2*(1/dx^2+1/dy^2);
        AA(kk,kk+nx) = 1/dy^2; % north
        AA(kk,kk-nx) = 1/dy^2; % south
        AA(kk,kk-1) = 1/dx^2; % west
        AA(kk,kk+1) = 1/dx^2; % east
        bb(kk) = -qqq/cond;
    end
end

% --- condizioni al contorno ---
% north (con north-west e north-west)
for jx = 1:nx
    kk = jx+(ny-1)*nx;
    AA(kk,kk) = 1;
    bb(kk) = Tnorth;
end

% south (con south-west e south-east)
for jx = 1:nx
    kk = jx;
    AA(kk,kk) = 1;
    bb(kk) = Tsouth;
end

% west (senza north-west e south-west)
for jy = 2:ny-1
    kk = 1+(jy-1)*nx;
    AA(kk,kk) = -cond*(dx/dy+dy/dx);
    AA(kk,kk+nx) = cond*dx/2/dy; % north
    AA(kk,kk-nx) = cond*dx/2/dy; % south
    AA(kk,kk+1) = cond*dy/dx; % east
    bb(kk) = -qqwest*dy;
end

% east (senza north-east e south-east)
for jy = 2:ny-1
    kk = jy*nx;
    AA(kk,kk) = -cond*(dx/dy+dy/dx)-hheast*dy;
    AA(kk,kk+nx) = cond*dx/2/dy; % north
    AA(kk,kk-nx) = cond*dx/2/dy; % south
    AA(kk,kk-1) = cond*dy/dx; % west
    bb(kk) = -hheast*Teast*dy;
end

%% risolvo e plotto
TT = AA\bb;
TTmat = reshape(TT,nx,ny);
surf(xmat,ymat,TTmat')
colorbar
xlabel('x [m]');
ylabel('y [m]');
zlabel('T [^oC]');
title('Stazionario 2D')

%% calcolo la potenza dispersa
qqlost = trapz(yvet,hheast*(TTmat(nx,:)-Teast));
fprintf('Calore perso: %.3f W/m\n',qqlost)