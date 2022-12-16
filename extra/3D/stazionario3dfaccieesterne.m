clc; clear; close all;
% il cubo in questione ha la faccia east che sottrae calore, la faccia west
% che scalda con convezione a 200 gradi, tutte le altre facce fissate con
% dirichlet a 0 gradi e una generazione di interna.
%% proprietà
qqq = 1e1;
cond = 2600; % esagerato

Lx = 1;
Ly = 1;
Lz = 1;

hhup = 100;
Tup = 20;

hhdown = 100;
Tdown = 20;

hhnorth = 100;
Tnorth = 20;

hhsouth = 100;
Tsouth = 20;

qqeast = -1e2;

hhwest = 100;
Twest = 100;

%% griglia spaziale
nn = 21; % nx | ny | nz
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
% spigoli up
for jx = 1:1
    for jy = 1:1
        kk = jx + (jy-1)*nn + (nn-1)*nn^2;
        
        AA(kk,kk) = 1;
        bb(kk) = Tup;
    end
end
for jx = 1:1
    for jy = 2:nn-1
        kk = jx + (jy-1)*nn + (nn-1)*nn^2;
        
        AA(kk,kk) = -dx/2*dy/dz*2-dz/2*dx/2/dy*2-dy*dz/2/dx*2-hhup/cond*dy/dz*2*dx/2-hhwest/cond*dy*dz/2/dx*2;
        AA(kk,kk-nn^2) = dx/2*dy/dz*2; % down
        AA(kk,kk+nn) = dx/2*dz/2/dy; % north
        AA(kk,kk-nn) = dx/2*dz/2/dy; % south
        AA(kk,kk+1) = dy*dz/2/dx*2; % east
        bb(kk) = -Tup*hhup/cond*dy/dz*2*dx/2-Twest*hhwest/cond*dy*dz/2/dx*2;
    end
end
for jx = 1:1
    for jy = nn:nn
        kk = jx + (jy-1)*nn + (nn-1)*nn^2;
        
        AA(kk,kk) = 1;
        bb(kk) = Tup;
    end
end
for jx = nn:nn
    for jy = 1:1
        kk = jx + (jy-1)*nn + (nn-1)*nn^2;
        
        AA(kk,kk) = 1;
        bb(kk) = Tup;
    end
end
for jx = nn:nn
    for jy = 2:nn-1
        kk = jx + (jy-1)*nn + (nn-1)*nn^2;
        
        AA(kk,kk) = -dx/2*dy/dz*2-dz/2*dx/2/dy*2-dy*dz/2/dx*2-hhup/cond*dy/dz*2*dx/2;
        AA(kk,kk-nn^2) = dx/2*dy/dz*2; % down
        AA(kk,kk+nn) = dx/2*dz/2/dy; % north
        AA(kk,kk-nn) = dx/2*dz/2/dy; % south
        AA(kk,kk-1) = dy*dz/2/dx*2; % west
        bb(kk) = -Tup*hhup/cond*dy/dz*2*dx/2-qqeast/cond*dy*dz/2/dx*2;
    end
end
for jx = nn:nn
    for jy = nn:nn
        kk = jx + (jy-1)*nn + (nn-1)*nn^2;
        
        AA(kk,kk) = 1;
        bb(kk) = Tup;
    end
end
for jx = 1:1
    for jy = 1:1
        kk = jx + (jy-1)*nn + (nn-1)*nn^2;
        
        AA(kk,kk) = 1;
        bb(kk) = Tup;
    end
end
for jx = 1:nn
    for jy = 1:1
        kk = jx + (jy-1)*nn + (nn-1)*nn^2;
        % sbagliato
        AA(kk,kk) = -dx/2*dy/dz*2-dz/2*dx/2/dy*2-dy*dz/2/dx*2-hhsouth/cond*dy/dz*2*dx/2-hhwest/cond*dy*dz/2/dx*2;
        AA(kk,kk-nn^2) = dx/2*dy/dz*2; % down
        AA(kk,kk+nn) = dx/2*dz/2/dy; % north
        AA(kk,kk-1) = dy*dz/2/dx*2; % west
        AA(kk,kk+1) = dy*dz/2/dx*2; % east
        bb(kk) = -Tup*hhup/cond*dy/2/dz*2*dx-Tsouth*hhsouth/cond*dy/2*dz/2/dx;
    end
end
for jx = 1:nn
    for jy = nn:nn
        kk = jx + (jy-1)*nn + (nn-1)*nn^2;
        
        AA(kk,kk) = 1;
        bb(kk) = Tup;
    end
end
for jx = 1:nn
    for jy = nn:nn
        kk = jx + (jy-1)*nn + (nn-1)*nn^2;
        
        AA(kk,kk) = 1;
        bb(kk) = Tup;
    end
end

% up (senza spigoli) - robin
for jx = 2:nn-1
    for jy = 2:nn-1
        kk = jx + (jy-1)*nn + (nn-1)*nn^2;
        
        AA(kk,kk) = -dx*dy/dz*2-dz*dx/dy-dy*dz/dx-hhup/cond*dy/dz*2*dx;
        AA(kk,kk-nn^2) = dx*dy/dz*2; % down
        AA(kk,kk+nn) = dx*dz/2/dy; % north
        AA(kk,kk-nn) = dx*dz/2/dy; % south
        AA(kk,kk+1) = dy*dz/2/dx; % east
        AA(kk,kk-1) = dy*dz/2/dx; % west
        bb(kk) = -Tup*hhup/cond*dy*dz/dx*2;
    end
end

% spigoli down
for jx = 1:1
    for jy = 1:nn
        kk = jx + (jy-1)*nn;
        
        AA(kk,kk) = 1;
        bb(kk) = Tdown;
    end
end
for jx = nn:nn
    for jy = 1:nn
        kk = jx + (jy-1)*nn;
        
        AA(kk,kk) = 1;
        bb(kk) = Tdown;
    end
end
for jx = 1:nn
    for jy = 1:1
        kk = jx + (jy-1)*nn;
        
        AA(kk,kk) = 1;
        bb(kk) = Tdown;
    end
end
for jx = 1:nn
    for jy = nn:nn
        kk = jx + (jy-1)*nn;
        
        AA(kk,kk) = 1;
        bb(kk) = Tdown;
    end
end

% down (senza spigoli) - robin
for jx = 2:nn-1
    for jy = 2:nn-1
        kk = jx + (jy-1)*nn;
        
        AA(kk,kk) = -dx*dy/dz*2-dz*dx/dy-dy*dz/dx-hhdown/cond*dy/dz*2*dx;
        AA(kk,kk+nn^2) = dx*dy/dz*2; % up
        AA(kk,kk+nn) = dx*dz/2/dy; % north
        AA(kk,kk-nn) = dx*dz/2/dy; % south
        AA(kk,kk+1) = dy*dz/2/dx; % east
        AA(kk,kk-1) = dy*dz/2/dx; % west
        bb(kk) = -Tdown*hhdown/cond*dy*dz/dx*2;
    end
end

% spigoli north
for jx = 1:1
    for jz = 2:nn-1
        kk = jx + (nn-1)*nn + (jz-1)*nn^2;
        
        AA(kk,kk) = 1;
        bb(kk) = Tnorth;
    end
end
for jx = nn:nn
    for jz = 2:nn-1
        kk = jx + (nn-1)*nn + (jz-1)*nn^2;
        
        AA(kk,kk) = 1;
        bb(kk) = Tnorth;
    end
end

% north (senza spigoli) - robin
for jx = 2:nn-1
    for jz = 2:nn-1
        kk = jx + (nn-1)*nn + (jz-1)*nn^2;
        
        AA(kk,kk) = -dx*dy/dz-dz*dx/dy*2-dy*dz/dx-hhnorth/cond/dy*dz*2*dx;
        AA(kk,kk+nn^2) = dx*dy/2/dz; % up
        AA(kk,kk-nn^2) = dx*dy/2/dz; % down
        AA(kk,kk-nn) = dx*dz/dy*2; % south
        AA(kk,kk+1) = dy*dz/2/dx; % east
        AA(kk,kk-1) = dy*dz/2/dx; % west
        bb(kk) = -Tnorth*hhnorth/cond/dy*dz*dx*2;
    end
end

% spigoli south
for jx = 1:1
    for jz = 2:nn-1
        kk = jx + (jz-1)*nn^2;
        
        AA(kk,kk) = 1;
        bb(kk) = Tnorth;
    end
end
for jx = nn:nn
    for jz = 2:nn-1
        kk = jx + (jz-1)*nn^2;
        
        AA(kk,kk) = 1;
        bb(kk) = Tnorth;
    end
end

% south (con spigoli west-east) - robin
for jx = 2:nn-1
    for jz = 2:nn-1
        kk = jx + (jz-1)*nn^2;
        
        AA(kk,kk) = -dx*dy/dz-dz*dx/dy*2-dy*dz/dx-hhsouth/cond/dy*dz*2*dx;
        AA(kk,kk+nn^2) = dx*dy/2/dz; % up
        AA(kk,kk-nn^2) = dx*dy/2/dz; % down
        AA(kk,kk+nn) = dx*dz/dy*2; % north
        AA(kk,kk+1) = dy*dz/2/dx; % east
        AA(kk,kk-1) = dy*dz/2/dx; % west
        bb(kk) = -Tsouth*hhsouth/cond/dy*dz*dx*2;
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
% Down
subplot(2,3,1)
surf(reshape(xmat(:,:,1),nn,nn),reshape(ymat(:,:,1),nn,nn),reshape(zmat(:,:,1),nn,nn),reshape(TTmat(:,:,1),nn,nn),'facecolor','interp')
title('Down')
colorbar
xlabel('x direction [m]')
ylabel('y direction [m]')
zlabel('z direction [m]')
grid on 
box on
set(gca,'fontsize',14)

% up
subplot(2,3,2)
surf(reshape(xmat(:,:,nn),nn,nn),reshape(ymat(:,:,nn),nn,nn),reshape(zmat(:,:,nn),nn,nn),reshape(TTmat(:,:,nn),nn,nn),'facecolor','interp')
title('Up')
colorbar
xlabel('x direction [m]')
ylabel('y direction [m]')
zlabel('z direction [m]')
grid on 
box on
set(gca,'fontsize',14)

% north
subplot(2,3,3)
surf(reshape(xmat(:,nn,:),nn,nn),reshape(ymat(:,nn,:),nn,nn),reshape(zmat(:,nn,:),nn,nn),reshape(TTmat(:,nn,:),nn,nn),'facecolor','interp')
title('North')
colorbar
xlabel('x direction [m]')
ylabel('y direction [m]')
zlabel('z direction [m]')
grid on 
box on
set(gca,'fontsize',14)

% south
subplot(2,3,4)
surf(reshape(xmat(:,1,:),nn,nn),reshape(ymat(:,1,:),nn,nn),reshape(zmat(:,1,:),nn,nn),reshape(TTmat(:,1,:),nn,nn),'facecolor','interp')
title('South')
colorbar
xlabel('x direction [m]')
ylabel('y direction [m]')
zlabel('z direction [m]')
grid on 
box on
set(gca,'fontsize',14)

% west
subplot(2,3,5)
surf(reshape(xmat(1,:,:),nn,nn),reshape(ymat(1,:,:),nn,nn),reshape(zmat(1,:,:),nn,nn),reshape(TTmat(1,:,:),nn,nn),'facecolor','interp')
title('West')
colorbar
xlabel('x direction [m]')
ylabel('y direction [m]')
zlabel('z direction [m]')
grid on 
box on
set(gca,'fontsize',14)

% east
subplot(2,3,6)
surf(reshape(xmat(nn,:,:),nn,nn),reshape(ymat(nn,:,:),nn,nn),reshape(zmat(nn,:,:),nn,nn),reshape(TTmat(nn,:,:),nn,nn),'facecolor','interp')
title('East')
colorbar
xlabel('x direction [m]')
ylabel('y direction [m]')
zlabel('z direction [m]')
grid on 
box on
set(gca,'fontsize',14)

set(gcf,'position',[200,100,2000,700])