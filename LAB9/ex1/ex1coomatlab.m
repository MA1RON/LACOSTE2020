clc; clear; close all;
%% dati
% topologia
Lx1 = .2;
Lx2 = .1;
Ly1 = .1;
Ly2 = .05;

% fisici
Text = 5;
Tint = 20;
cond = 10;

%% griglie
nx1 = 21;
nx2 = 11;
ny1 = 21;
ny2 = 11;

ntot = nx1*ny1-(nx1-nx2)*(ny1-ny2);
xvec = linspace(0,Lx1,nx1);
yvec = linspace(0,Ly1,ny1);
dx = xvec(2)-xvec(1);
dy = yvec(2)-yvec(1);
[xmat, ymat] = meshgrid(xvec,yvec);

%% creo il sistema
AA = sparse([],[],[],ntot,ntot,5*ntot);
bb = zeros(ntot,1);

% --- center ---
% center west block
for jx = 2:nx1-1
    for jy = 2:ny2-1
        kk = jx+(jy-1)*nx1;
        AA(kk,kk) = -2*cond*(1/dx^2+1/dy^2);
        AA(kk,kk-1) = cond/dx^2;
        AA(kk,kk+1) = cond/dx^2;
        AA(kk,kk-nx1) = cond/dy^2;
        AA(kk,kk+nx1) = cond/dy^2;
    end
end

% center separation
jy = nx2;
for jx = 2:nx2-1
    kk = jx + nx1*(jy-1);
    AA(kk) = -2*cond*(1/dx^2+1/dy^2);
    AA(kk,kk-1) = cond/dx^2; % north
    AA(kk,kk+1) = cond/dx^2; % south
    AA(kk,kk-nx1) = cond/dy^2; % west
    AA(kk,kk+nx1) = cond/dy^2; % east
end

% center east block
for jx = 2:nx2-1
    for jy = ny2+1:ny1-1
        kk = jx+(jy-ny2-1)*nx2 + nx1*ny2;
        AA(kk) = -2*cond*(1/dx^2+1/dy^2);
        AA(kk,kk-1) = cond/dx^2; % north
        AA(kk,kk+1) = cond/dx^2; % south
        if jy == ny2+1
            AA(kk,kk-nx1) = cond/dy^2; % west
        else
            AA(kk,kk-nx2) = cond/dy^2; % west
        end
        AA(kk,kk+nx2) = cond/dy^2; % east
    end
end

% --- condizioni al contorno ---
% north
jx = 1;
for jy = 1:ny1
    if jy<=nx2
        kk = jx+(jy-1)*nx1;
    else
        kk = jx+(jy-1-ny2)*nx2+nx1*ny2;
    end
    AA(kk,kk) = 1;
    bb(kk) = Text;
end

% south esterno
jx = nx1;
for jy = 2:nx2-1
    kk = jx+(jy-1)*nx1;
    %{
    AA(kk,kk) = -2*cond*(1/dx^2+1/dy^2);
    AA(kk,kk-1) = 2*cond/dx^2; % north
    AA(kk,kk-nx1) = 1*cond/dy^2; % west
    AA(kk,kk+nx1) = 1*cond/dy^2; % east
    %}
    AA(kk,kk) = 1;
    bb(kk) = Text;
end

% south interno
jx = nx2;
for jy = ny2+1:ny1
    kk = nx1*ny2+jx+(jy-1-ny2)*nx2;
    AA(kk,kk) = 1;
    bb(kk) = Tint;
end

% west
for jx = 2:nx1
    kk = jx;
    AA(kk,kk) = 1;
    bb(kk) = Text;
end

% east interno
jy = ny2;
for jx = nx2:nx1
    kk = jx + (jy-1)*nx1;
    AA(kk,kk) = 1;
    bb(kk) = Tint;
end

% east esterno
jy = nx1;
for jx = 2:ny2-1
    kk = jx+(jy-1-ny2)*nx2+nx1*ny2;
    %{
    AA(kk,kk) = -2*cond*(1/dx^2+1/dy^2);
    AA(kk,kk-1) = cond/dx^2; % north
    AA(kk,kk+1) = cond/dx^2; % south
    AA(kk,kk-nx2) = 2*cond/dy^2; % west
    %}
    AA(kk,kk) = 1;
    bb(kk) = Text;
end

% --- risolvo ---

TT = AA\bb;

TTmat = zeros(nx1,ny1);
for jx = 1:nx1
    for jy = 1:ny1
        if jx>nx2 && jy>ny2
            TTmat(jx,jy) = nan;
        elseif jy<=ny2
            kk = jx+(jy-1)*nx1;
            TTmat(jx,jy) = TT(kk);
        else
            kk = jx+(jy-1-ny2)*nx2+nx1*ny2;
            TTmat(jx,jy) = TT(kk);
        end
    end
end

%% surfo 
surf(xmat, ymat, TTmat, 'facecolor', 'interp')