clc; clear; close all;
%% proprietà
qqq = 5e3;
cond = 1e4; % esagerato

Lx = 1;
Ly = 2;
Lz = 1;

Tup = 0;
Tdown = 0;
Tnorth = 0;
Tsouth = 0;
qqeast = -1e4; % refrigera
hhwest = 1e4;
Twest = 200;

%% griglia spaziale
nn = 11; % nx | ny | nz
ntot = nn^3;

xvet = linspace(0,Lx,nn);
yvet = linspace(0,Ly,nn);
zvet = linspace(0,Lz,nn);

dx = xvet(2)-xvet(1);
dy = yvet(2)-yvet(1);
dz = zvet(2)-zvet(1);

[xmat,ymat,zmat] = meshgrid(xvet,yvet,zvet);

%% creo il problema
AA = sparse([],[],[],ntot,ntot,7*ntot);
bb = zeros(ntot,1);

% --- center block ---
for jx = 2:nn-1
    for jy = 2:nn-1
        for jz = 2:nn-1
            kk = jx + (jy-1)*nn + (jz-1)*nn^2;
            
            AA(kk,kk) = -2*(dx*dy/dz + dy*dz/dx + dz*dx/dy);
            AA(kk,kk+nn^2) = dx*dy/dz; % up
            AA(kk,kk-nn^2) = dx*dy/dz; % down
            AA(kk,kk+nn) = dz*dx/dy; % north
            AA(kk,kk-nn) = dz*dx/dy; % south
            AA(kk,kk+1) = dy*dz/dx; % east
            AA(kk,kk-1) = dy*dz/dx; % west
            bb(kk) = - qqq/cond;
        end
    end
end

% --- condizioni al contorno ---
% up (con spigoli) - dirichlet
for jx = 1:nn
    for jy = 1:nn
        kk = jx + (jy-1)*nn + (nn-1)*nn^2;
        
        AA(kk,kk) = 1;
        bb(kk) = Tup;
    end
end

% down (con spigoli) - dirichlet
for jx = 1:nn
    for jy = 1:nn
        kk = jx + (jy-1)*nn;
        
        AA(kk,kk) = 1;
        bb(kk) = Tdown;
    end
end

% north (con spigoli west-east) - dirichlet
for jx = 1:nn
    for jz = 2:nn-1
        kk = jx + (nn-1)*nn + (jz-1)*nn^2;
        
        AA(kk,kk) = 1;
        bb(kk) = Tnorth;
    end
end

% south (con spigoli west-east) - dirichlet
for jx = 1:nn
    for jz = 2:nn-1
        kk = jx + (jz-1)*nn^2;
        
        AA(kk,kk) = 1;
        bb(kk) = Tsouth;
    end
end

% east (senza spigoli) - neumann non omogeneo
for jy = 2:nn-1
    for jz = 2:nn-1
        kk = nn + (jy-1)*nn + (jz-1)*nn^2;
        
        AA(kk,kk) = -(dx*dy/dz+dz*dx/dy)-dy*dz/dx*2;
        AA(kk,kk+nn^2) = dx/2*dy/dz; % up
        AA(kk,kk-nn^2) = dx/2*dy/dz; % down
        AA(kk,kk+nn) = dx/2*dz/dy; % north
        AA(kk,kk-nn) = dx/2*dz/dy; % south
        AA(kk,kk-1) = dy*dz/dx*2; % west
        bb(kk) = - qqeast/cond;
    end
end

% west (senza spigoli) - robin
for jy = 2:nn-1
    for jz = 2:nn-1
        kk = 1 + (jy-1)*nn + (jz-1)*nn^2;
        
        AA(kk,kk) = -(dx*dy/dz+dz*dx/dy)-dy*dz/dx*2-hhwest/cond*dy*dz/dx*2;
        AA(kk,kk+nn^2) = dx/2*dy/dz; % up
        AA(kk,kk-nn^2) = dx/2*dy/dz; % down
        AA(kk,kk+nn) = dx/2*dz/dy; % north
        AA(kk,kk-nn) = dx/2*dz/dy; % south
        AA(kk,kk+1) = dy*dz/dx*2; % east
        bb(kk) = -Twest*hhwest/cond*dy*dz/dx*2;
    end
end

%% risolvo il sistema
TT = AA\bb;

TTmat = reshape(TT,nn,nn,nn);

%% post production
for jx = 1:nn
    for jy = 1:nn
        for jz = 1:nn
            kk = jx + (jy-1)*nn + (jz-1)*nn^2;
            scatter3(jx,jy,jz,20,TTmat(kk));
            hold on
        end
    end
end
colorbar
title('Stazionario 3D')
xlabel('x direction [m]')
ylabel('y direction [m]')
zlabel('z direction [m]')
grid on 
box on
set(gca,'fontsize',16)