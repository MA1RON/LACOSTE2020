clc; clear; close all;
%% dati
LL = 1;
Tnorth = 500;
Tsouth = 500;
Twest = 500;
Teast = 300;
hheast = 10;
cond = .72;

%% griglia spaziale
nx = 51;
ny = 51;
ntot = nx*ny;
xvet = linspace(0,LL,nx);
yvet = linspace(0,LL,ny);
dx = xvet(2)-xvet(1);
dy = yvet(2)-yvet(1);
[xmat,ymat] = meshgrid(xvet,yvet);

%% creo il problema
AA = sparse([],[],[],ntot,ntot,5*ntot);
bb = zeros(ntot,1);

% --- centro ---
for jx = 2:nx-1
    for jy = 2:ny-1
        kk = jx+(jy-1)*nx;
        AA(kk,kk) = -2*(1/dx^2+1/dy^2);
        AA(kk,kk+nx) = 1/dy^2; % north
        AA(kk,kk-nx) = 1/dy^2; % south
        AA(kk,kk+1) = 1/dx^2; % west
        AA(kk,kk-1) = 1/dx^2; % east
    end
end

% --- condizioni al contorno ---
% north (con north-west e north-east) - dirichlet
for jx = 1:nx
    kk = jx+(ny-1)*nx;
    AA(kk,kk) = 1;
    bb(kk) = Tnorth;
end

% south (con south-west e south-east) - dirichlet
for jx = 1:nx
    kk = jx;
    AA(kk,kk) = 1;
    bb(kk) = Tsouth;
end

% west (senza north-west e south-west) - dirichlet
for jy = 2:ny-1
    kk = 1+(jy-1)*nx;
    AA(kk,kk) = 1;
    bb(kk) = Twest;
end

% east (senza north-east e south-east) - robin
for jy = 2:ny-1
    kk = jy*nx;
    AA(kk,kk) = -cond*(dx/dy+dy/dx)-hheast*dy;
    AA(kk,kk+nx) = cond*dx/2/dy; % north
    AA(kk,kk-nx) = cond*dx/2/dy; % south
    AA(kk,kk-1) = cond*dy/dx; % west
    bb(kk) = -hheast*Teast*dy;
end

%% --- risolvo ---
TT = AA\bb;
TTmat = reshape(TT,nx,ny);

%% calcolo della potenza dispersa ad east
% - da me -
qqeast = (TTmat(nx-1,1)-TTmat(nx,1))*cond*dy/(dx/2);
for jy = 2:ny-1
    qqeast = qqeast + (TTmat(nx-1,jy)-TTmat(nx,jy))*cond*dy/dx;
end
qqeast = qqeast + (TTmat(nx-1,ny)-TTmat(nx,ny))*cond*dy/(dx/2);
disp(['Calore conduttivo (approssimazione) = ',num2str(qqeast),' [W/m]']);

% - da prof -
qlin=trapz(yvet,hheast*(TT(nx:nx:end)-Teast));
disp(['Calore convettivo (corretto) = ',num2str(qlin),' [W/m]']);

%% post production - temperatura stazionaria 2D
surf(xmat,ymat,TTmat-273.15,'facecolor','interp')
xlabel('Direzione x [m]')
ylabel('Direzione y [m]')
zlabel('Temperatura [^oC]')
set(gca,'fontsize',18)
title('Stazionario 2D')