clc; clear; close all;
%% dati
% geometrici
Lx1 = .1;
Lx2 = .2;
Ly1 = .05;
Ly2 = .1;

% fisici
Tint = 20;
Text = 5;
cond = 10;

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
        % bb(kk) = 0 non ho generazione interna
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
        % bb(kk) = 0 non ho generazione interna
    end
end

% --- condizioni al contorno ---
% north (con north-west and north-east)
for jx = 1:nx2
    kk = nx1*(ny1-1)+nx2*(ny2-ny1)+jx;
    AA(kk,kk) = 1;
    bb(kk) = Text;
end

% south esterno
for jx = 2:nx1-1
    kk = jx;
    AA(kk,kk) = -2*(1/dx^2+1/dy^2);
    AA(kk,kk+nx1) = 2/dy^2; % north
    AA(kk,kk-1) = 1/dx^2; % west
    AA(kk,kk+1) = 1/dx^2; % east
    % bb(kk) = 0 neumann omogeneo
end

% south interno (con nodo centrale and east-center)
for jx = nx1:nx2
    kk = nx1*(ny1-1)+jx;
    AA(kk,kk) = 1;
    bb(kk) = Tint;
end

% west (con south-west)
for jy = 1:ny2-1
    if jy<=ny1
        kk = (jy-1)*nx1+1;
    else
        kk = (ny1-1)*nx1+(jy-ny1)*nx2+1;
    end
    AA(kk,kk) = 1;
    bb(kk) = Text;
end

% east esterno
for jy = ny1+1:ny2-1
    kk = nx1*(ny1-1)+(jy-ny1+1)*nx2;
    AA(kk,kk) = -2*(1/dx^2+1/dy^2);
    AA(kk,kk+nx2) = 1/dy^2; % north
    AA(kk,kk-nx2) = 1/dy^2; % south
    AA(kk,kk-1) = 2/dx^2; % west
    % bb(kk) = 0 neumann omogeneo per simmetria
end

% east interno (con south-center)
for jy = 1:ny1-1
    kk = nx1+(jy-1)*nx1;
    AA(kk,kk) = 1;
    bb(kk) = Tint;
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