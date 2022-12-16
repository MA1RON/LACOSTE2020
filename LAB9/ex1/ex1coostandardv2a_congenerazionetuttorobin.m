clc; clear; close all;
%% dati
% geometrici
Lx1 = .1;
Lx2 = .2;
Ly1 = .05;
Ly2 = .1;

% fisici
Tint = 20;
hhint = 50;
Text = 5;
hhext = 1000; % vento molto forte
cond = 10;
qqq = 1e4; % solo in mezzo al muro per tubi passanti

%% griglia
nx2 = 51;
nx1 = ceil(nx2/2);
ny2 = 21;
ny1 = ceil(ny2/2);
ntot = nx1*(ny1-1)+nx2*(ny2-ny1+1);

xvet = linspace(0,Lx2,nx2);
yvet = linspace(0,Ly2,ny2);
dx = xvet(2)-xvet(1);
dy = yvet(2)-yvet(1);

[xmat,ymat] = meshgrid(xvet,yvet);

%% creo il problema
AA = sparse([],[],[],ntot,ntot,5*ntot);
bb = zeros(ntot,1);

% --- center ---
% center block south
for jx = 2:nx1-1
    for jy = 2:ny1
        kk = jx+(jy-1)*nx1;
        AA(kk,kk) = -2*cond*(1/dx^2+1/dy^2);
        if not(jy==ny1)
            AA(kk,kk+nx1) = cond/dy^2; % north
        else
            AA(kk,kk+nx2) = cond/dy^2; % north
        end
        AA(kk,kk-nx1) = cond/dy^2; % south
        AA(kk,kk-1) = cond/dx^2; % west
        AA(kk,kk+1) = cond/dx^2; % east
        % generazione solo al centro
        if jx==round(nx1/2) 
            bb(kk) = -qqq/cond/dx;
        end
    end
end

% center block north
for jx = 2:nx2-1
    for jy = ny1+1:ny2-1
        kk = jx+(jy-ny1-1)*nx2+nx1*ny1+nx2-nx1;
        AA(kk,kk) = -2*cond*(1/dx^2+1/dy^2);
        AA(kk,kk+nx2) = cond/dy^2; % north
        AA(kk,kk-nx2) = cond/dy^2; % south
        AA(kk,kk-1) = cond/dx^2; % west
        AA(kk,kk+1) = cond/dx^2; % east
        % generazione solo al centro
        if jy==round(ny1/2*3) && jx>=round(nx1/2) 
            bb(kk) = -qqq/cond/dy;
        elseif jx==round(nx1/2) && jy<=round(ny1/2*3)
            bb(kk) = -qqq/cond/dx;
        end
    end
end

% --- condizioni al contorno ---
% south esterno - neumann
for jx = 2:nx1-1
    kk = jx;
    AA(kk,kk) = -2*(1/dx^2+1/dy^2);
    AA(kk,kk+nx1) = 2/dy^2; % north
    AA(kk,kk-1) = 1/dx^2; % west
    AA(kk,kk+1) = 1/dx^2; % east
    % bb(kk) = 0 neumann omogeneo
end

% south west - robin da west e nemaunn da south
for jy = 1:1
    kk = 1;
    AA(kk,kk) = -(dx/dy+dy/dx)*cond/2-dy/2*hhext;
    AA(kk,kk+nx1) = dx/2*cond/dy; % north
    AA(kk,kk+1) = dy/2*cond/dx; % east
    bb(kk) = -dy/2*hhext*Text;
end

% west - robin
for jy = 2:ny2-1
    if jy <= ny1
        kk = (jy-1)*nx1+1;
        ks = kk - nx1;
    else
        kk = (ny1-1)*nx1+(jy-ny1)*nx2+1;
        ks = kk - nx2;
    end
    if jy < ny1
        kn = kk + nx1;
    else
        kn = kk + nx2;
    end
    AA(kk,kk) = -(dx/dy+dy/dx)*cond-dy*hhext;
    AA(kk,kn) = dx/2*cond/dy; % north
    AA(kk,ks) = dx/2*cond/dy; % south
    AA(kk,kk+1) = dy*cond/dx; % east
    bb(kk) = -dy*hhext*Text;
end

% CONDIZIONE INCRIMINATA
% north west - robin da north e robin da west
for jx = 1:1
    kk = nx1*(ny1-1)+nx2*(ny2-ny1)+jx;
    AA(kk,kk) = -(dx/dy+dy/dx)*cond/2-(dx+dy)/2*hhext;
    AA(kk,kk-nx2) = dx/2*cond/dy; % south
    AA(kk,kk+1) = dy/2*cond/dx; % east
    bb(kk) = -(dy+dx)/2*hhext*Text;
end

% north - robin
for jx = 1:nx2-1
    kk = nx1*(ny1-1)+nx2*(ny2-ny1)+jx;
    AA(kk,kk) = -(dx/dy+dy/dx)*cond-dx*hhext;
    AA(kk,kk-nx2) = dx*cond/dy; % south
    AA(kk,kk-1) = dy/2*cond/dx; % west
    AA(kk,kk+1) = dy/2*cond/dx; % east
    bb(kk) = -dx*hhext*Text;
end

% north-east - robin da north e neumann da east
for jx = nx2:nx2
    kk = nx1*(ny1-1)+nx2*(ny2-ny1)+jx;
    AA(kk,kk) = -(dx/dy+dy/dx)*cond/2-dx/2*hhext;
    AA(kk,kk-nx2) = dx/2*cond/dy; % south
    AA(kk,kk-1) = dy/2*cond/dx; % west
    bb(kk) = -dx/2*hhext*Text;
end

% east esterno - neumann omogeneo
for jy = ny1+1:ny2-1
    kk = nx1*(ny1-1)+(jy-ny1+1)*nx2;
    AA(kk,kk) = -2*(1/dx^2+1/dy^2);
    AA(kk,kk+nx2) = 1/dy^2; % north
    AA(kk,kk-nx2) = 1/dy^2; % south
    AA(kk,kk-1) = 2/dx^2; % west
    % bb(kk) = 0 neumann omogeneo per simmetria
end

% east-center
for jx = nx2:nx2
    kk = nx1*(ny1-1)+jx;
    AA(kk,kk) = -(dx/dy+dy/dx)/2*cond-dy*hhint/2;
    AA(kk,kk+nx2) = dx/2*cond/dy; % north
    AA(kk,kk-1) = dy/2*cond/dx; % west
    bb(kk) = -dy*hhint/2*Tint;
end

% south interno - robin
for jx = nx1+1:nx2-1
    kk = nx1*(ny1-1)+jx;
    AA(kk,kk) = -dx*cond/dy-dy*cond/dx-dx*hhint;
    AA(kk,kk+nx2) = dx*cond/dy; % north
    AA(kk,kk-1) = dy*cond/2/dx; % west
    AA(kk,kk+1) = dy*cond/2/dx; % east
    bb(kk) = -dx*hhint*Tint;
end

% nodo centrale - robin
for jx = nx1:nx1
    kk = nx1*(ny1-1)+jx;
    AA(kk,kk) = -3/2*dx*cond/dy-3/2*dy*cond/dx-dx/2*hhint-dy*hhint/2;
    AA(kk,kk+nx2) = dx*cond/dy; % north
    AA(kk,kk-nx1) = dx*cond/2/dy; % south
    AA(kk,kk-1) = dy*cond/dx; % west
    AA(kk,kk+1) = dy*cond/dx/2; % east
    bb(kk) = -(dx+dy)/2*hhint*Tint;
end

% east interno - robin da east e neumann da south
for jy = 2:ny1-1
    kk = nx1+(jy-1)*nx1;
    AA(kk,kk) = -dx*cond/dy-dy*cond/dx-dy*hhint;
    AA(kk,kk+nx1) = dx*cond/dy/2; % north
    AA(kk,kk-nx1) = dx*cond/dy/2; % south
    AA(kk,kk-1) = dy*cond/dx; % west
    bb(kk) = -dy*hhint*Tint;
end

% south center - robin
for jy = 1:1
    kk = nx1+(jy-1)*nx1;
    AA(kk,kk) = -dx/2*cond/dy-dy*cond/2/dx-dy*hhint/2;
    AA(kk,kk+nx1) = dx/2*cond/dy; % north
    AA(kk,kk-1) = dy/2*cond/dx; % west
    bb(kk) = -dy*hhint/2*Tint;
end

% --- risolvo ---
TT = AA\bb;

%% reshape
TTmat = zeros(nx2,ny2);
for jx = 1:nx2
    for jy = 1:ny2
        if jx>nx1 && jy<ny1
            TTmat(jx,jy) = nan;
        else
            if jy<ny1
                kk = jx+(jy-1)*nx1;
            else
                kk = jx+(jy-ny1)*nx2+nx1*(ny1-1);
            end
            TTmat(jx,jy) = TT(kk);
        end
    end
end

%% post production
surf(xmat,ymat,TTmat','facecolor','interp');
xlabel('Direzione x [m]')
ylabel('Direzione y [m]')
zlabel('Temperatura [^oC]')
set(gca,'fontsize',18)
title('Stazionario 2D')