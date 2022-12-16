clc; clear; close all;
%% dati
Ly = .4;
Lx = .2;
cond = 400;
Tinf = 25;
hh = 3e4;
Tnorth = 75;
Lchip = 10e-3;
Achip = Lchip^2;

%% griglie spaziali
nx2 = 31;
nx1 = round(nx2/2);
ny2 = 31;
ny1 = round(ny2/2);
ntot = nx2*ny2;

xvet = linspace(0,Lx,nx2);
yvet = linspace(0,Ly,ny2);
dx = xvet(2)-xvet(1);
dy = yvet(2)-yvet(1);

[xmat,ymat] = meshgrid(xvet,yvet);

%% creo il sistema
AA = sparse([],[],[],ntot,ntot,5*ntot);
bb = zeros(ntot,1);

% --- centro ---
% south block
for jx = 2:nx1-1
    for jy = 2:ny1 % compreso il bordo tra i due blocchi
        kk = jx+(jy-1)*nx1;
        
        AA(kk,kk) = -2*(dx/dy+dy/dx);
        if jy == ny1
            AA(kk,kk+nx2) = dx/dy; % north
        else
            AA(kk,kk+nx1) = dx/dy; % north
        end
        AA(kk,kk-nx1) = dx/dy; % south
        AA(kk,kk-1) = dy/dx; % west
        AA(kk,kk+1) = dy/dx; % east
    end
end

% north block
for jx = 2:nx2-1
    for jy = ny1+1:ny2-1
        kk = nx1*(ny1-1)+jx+(jy-ny1)*nx2;
        
        AA(kk,kk) = -2*(dx/dy+dy/dx);
        AA(kk,kk+nx2) = dx/dy; % north
        AA(kk,kk-nx2) = dx/dy; % south
        AA(kk,kk-1) = dy/dx; % west
        AA(kk,kk+1) = dy/dx; % east
    end
end

% --- condizioni al contorno ---
% north (con north-west e north-east)
for jx = 1:nx2
    kk = jx+(ny2-ny1)*nx2+(ny1-1)*nx1;
    
    AA(kk,kk) = 1;
    bb(kk) = Tnorth;
end

% east esterno (senza nort-east e center-east)
for jy = ny1+1:ny2-1
    kk = (ny1-1)*nx1+nx2*(jy-ny1+1);
    
    AA(kk,kk) = -(dx/dy+dy/dx);
    AA(kk,kk+nx2) = dx/2/dy; % north
    AA(kk,kk-nx2) = dx/2/dy; % south
    AA(kk,kk-1) = dy/dx; % west
end

% center-east (neumann da east e robin da south)
for jx = nx2:nx2
    kk = jx+(ny1-1)*nx1;
    
    AA(kk,kk) = -(dx/dy+dy/dx+hh/cond*dx)/2;
    AA(kk,kk+nx2) = dx/2/dy; % north
    AA(kk,kk-1) = dy/2/dx; % west
    bb(kk) = -hh/cond*dx/2*Tinf;
end

% south interno (senza center-east e center-center)
for jx = nx1+1:nx2-1
    kk = jx+(ny1-1)*nx1;
    
    AA(kk,kk) = -(dx/dy+dy/dx)-hh/cond*dx;
    AA(kk,kk+nx2) = dx/dy; % north
    AA(kk,kk-1) = dy/2/dx; % west
    AA(kk,kk+1) = dy/2/dx; % east
    bb(kk) = -hh/cond*dx*Tinf;
end

% center-center
for jx = nx1:nx1
    kk = jx*ny1;
    
    AA(kk,kk) = -dx/dy-dy/dx-dx/2/dy-dy/2/dx-hh/cond*(dx+dy)/2;
    AA(kk,kk+nx2) = dx/dy; % north
    AA(kk,kk-nx2) = dx/2/dy; % south
    AA(kk,kk-1) = dy/dx; % west
    AA(kk,kk+1) = dy/2/dx; % east
    bb(kk) = -hh/cond/2*Tinf*(dx+dy);
end

% east interno (senza center-center e south-center)
for jy = 2:ny1-1
    kk = jy*nx1;
    
    AA(kk,kk) = -(dx/dy+dy/dx)-hh/cond*dy;
    AA(kk,kk+nx1) = dx/2/dy; % north
    AA(kk,kk-nx1) = dx/2/dy; % south
    AA(kk,kk-1) = dy/dx; % west
    bb(kk) = -hh/cond*dy*Tinf;
end

% south-center (robin da east e neumann da south)
for jx = nx1:nx1
    kk = jx;
    
    AA(kk,kk) = -(dx/dy+dy/dx)/2-hh/cond*dy/2;
    AA(kk,kk+nx1) = dx/2/dy; % north
    AA(kk,kk-1) = dy/2/dx; % west
    bb(kk) = -hh/cond*Tinf*dy/2;
end

% south esterno (senza south-center e south-west)
for jx = 2:nx1-1
    kk = jx;
    
    AA(kk,kk) = -(dx/dy+dy/dx);
    AA(kk,kk+nx1) = dx/dy; % north
    AA(kk,kk-1) = dy/2/dx; % west
    AA(kk,kk+1) = dy/2/dx; % east
end

% south-west (neumann da east e da south)
for jx = 1:1
    kk = jx;
    
    AA(kk,kk) = -(dx/dy+dy/dx)/2;
    AA(kk,kk+nx1) = dx/2/dy; % north
    AA(kk,kk+1) = dy/2/dx; % east
end

% west (senza south-west e north-west)
for jy = 2:ny2-1
    if jy <= ny1
        kk = 1+(jy-1)*nx1;
    else
        kk = 1+(jy-ny1)*nx2+(ny1-1)*nx1;
    end
    
    AA(kk,kk) = -(dx/dy+dy/dx);
    if jy < ny1
        AA(kk,kk+nx1) = dx/2/dy; % north
    else
        AA(kk,kk+nx2) = dx/2/dy; % north
    end
    if jy <= ny1
        AA(kk,kk-nx1) = dx/2/dy; % south
    else
        AA(kk,kk-nx2) = dx/2/dy; % south
    end
    AA(kk,kk+1) = dy/dx; % east
end

% --- risolvo il sistema ---
TT = AA\bb;

%% sistemo il dominio
TTmat = zeros(nx2,ny2);
for jx = 1:nx2
    for jy = 1:ny2
        if jx>nx1 && jy<ny1
            TTmat(jx,jy) = nan;
        else
            if jy<=ny1
                kk = jx+(jy-1)*nx1;
            else
                kk = nx1*(ny1-1)+jx+(jy-ny1)*nx2;
            end
            TTmat(jx,jy) = TT(kk);
        end
    end
end

%% calcolo massima dissipazione di calore
% dissipazione south interna
nchipperdom = Lx/Lchip;
fprintf('Ci sono %i chip in ogni dominio considerato.\n',nchipperdom)
qconvsouth = trapz(xvet(nx1:nx2),hh*(Tinf-TTmat(nx1:nx2,ny1)));
qconveast = trapz(yvet(1:ny1-1),hh*(Tinf-TTmat(nx1,1:ny1-1)));
qconv = qconvsouth+qconveast;
fprintf('Calore disperso sul bordo orizzonatale: %.3f W.\n',-Lchip*qconvsouth/nchipperdom)
fprintf('Calore disperso sul bordo verticale: %.3f W.\n',-Lchip*qconveast/nchipperdom)
fprintf('Calore disperso totale: %.3f W.\n',-Lchip*qconv/nchipperdom)

%% post production
surf(xmat,ymat,TTmat','facecolor','interp')
colorbar