clc; clear; close all;
%% dati
% geometrici
Lx1 = .1;
Lx2 = .2;
Ly1 = .05;
Ly2 = .1;

% fisici
hh = 5;
Tint = 20;
Text = 5;
cond = 10;

%% griglia
nx2 = 101;
nx1 = ceil(nx2/2);
ny2 = 101;
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

% center-center
AA(nx1*(ny1-1)+nx1,nx1*(ny1-1)+nx1) = -(3/dx^2+3/dy^2+hh/cond*(1/dy+1/dx));
AA(nx1*(ny1-1)+nx1,nx1*(ny1-1)+nx1+nx2) = 2/dy^2; % north
AA(nx1*(ny1-1)+nx1,nx1*(ny1-1)+nx1-nx1) = 1/dy^2; % south
AA(nx1*(ny1-1)+nx1,nx1*(ny1-1)+nx1-1) = 2/dx^2; % west
AA(nx1*(ny1-1)+nx1,nx1*(ny1-1)+nx1+1) = 1/dx^2; % east
bb(nx1*(ny1-1)+nx1) = -Tint*hh/cond*(1/dx+1/dy);

% south interno
for jx = nx1+1:nx2-1
    kk = nx1*(ny1-1)+jx;
    AA(kk,kk) = -2*(1/dx^2+1/dy^2+hh/cond/dy);
    AA(kk,kk+nx2) = 2/dy^2; % north
    AA(kk,kk-1) = 1/dx^2; % west
    AA(kk,kk+1) = 1/dx^2; % east
    bb(kk) = -2*Tint*hh/cond/dy;
end

% east-center
AA(nx1*(ny1-1)+nx2,nx1*(ny1-1)+nx2) = -(1/dx^2+1/dy^2+hh/cond*(1/dy+1/dx));
AA(nx1*(ny1-1)+nx2,nx1*(ny1-1)+nx2+nx2) = 1/dy^2; % north
AA(nx1*(ny1-1)+nx2,nx1*(ny1-1)+nx2-1) = 1/dx^2; % west
bb(nx1*(ny1-1)+nx2) = -hh*Tint/cond*(1/dx+1/dy);

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

% south center
AA(nx1,nx1) = -(1/dx^2+1/dy^2+hh/cond*(1/dx+1/dy));
AA(nx1,nx1+nx1) = 1/dy^2; % north
AA(nx1,nx1-1) = 1/dx^2; % west
bb(nx1) = -hh*Tint/cond*(1/dy+1/dx);

% east interno
for jy = 2:ny1-1
    kk = nx1+(jy-1)*nx1;
    AA(kk,kk) = -2*(1/dx^2+1/dy^2+hh/cond/dx);
    AA(kk,kk+nx1) = 1/dy^2; % north
    AA(kk,kk-nx1) = 1/dy^2; % south
    AA(kk,kk-1) = 2/dx^2; % west
    bb(kk) = -2*hh*Tint/cond/dx;
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

%% calore disperso
% north
qnorth = (TTmat(1,ny2)*(1/dx^2+1/dy^2) - TTmat(2,ny2)/dx^2 - TTmat(1,ny2-1)/dy^2)*cond/(1/dx+1/dy);
jy = ny2;
for jx = 2:nx2-1
    qnorth = qnorth + 2*TTmat(jx,jy)*(1/dx^2+1/dy^2)*cond*dy/2; % center
    qnorth = qnorth - 2*TTmat(jx,jy-1)/dy^2*cond*dy/2; % south
    qnorth = qnorth - TTmat(jx-1,jy)/dx^2*cond*dy/2; % west
    qnorth = qnorth - TTmat(jx+1,jy)/dx^2*cond*dy/2; % east
end
qnorth = qnorth + (TTmat(nx2,ny2)*(1/dx^2+1/dy^2) - TTmat(nx2-1,ny2)/dx^2 - TTmat(nx2,ny2-1)/dy^2)*cond/(1/dx+1/dy);
fprintf('Il calore disperso a nord vale: %.3f W.\n', - qnorth)

% west
qwest = (TTmat(1,1)*(1/dx^2+1/dy^2) - TTmat(2,1)/dx^2 - TTmat(1,2)/dy^2)*cond/(1/dx+1/dy);
jx = 1;
for jy = 2:ny2-1
    qwest = qwest + 2*TTmat(jx,jy)*(1/dx^2+1/dy^2)*cond*dx/2; % center
    qwest = qwest - TTmat(jx,jy+1)/dy^2*cond*dx/2; % north
    qwest = qwest - TTmat(jx,jy-1)/dy^2*cond*dx/2; % south
    qwest = qwest - 2*TTmat(jx+1,jy)/dx^2*cond*dx/2; % east
end
qwest = qwest + (TTmat(1,ny2)*(1/dx^2+1/dy^2) - TTmat(2,ny2)/dx^2 - TTmat(1,ny2-1)/dy^2)*cond/(1/dx+1/dy);
fprintf('Il calore disperso a nord vale: %.3f W.\n', - qwest)


%% post production
figure
surf(xmat,ymat,TTmat','facecolor','interp');
xlabel('Direzione x [m]')
ylabel('Direzione y [m]')
zlabel('Temperatura [^oC]')
set(gca,'fontsize',18)
title('Stazionario 2D')