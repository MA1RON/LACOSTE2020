clc; clear; close all;
%% dati
Lx = 1;
Ly = 1;
cond = 1;

Ts = 10;
Tn = 20;
Te = 30;
Tw = 40;

%% griglie spaziali
nx = 100;
ny = 100;
ntot = nx*ny;
xvec = linspace(0,Lx,nx);
yvec = linspace(0,Ly,ny);
dx = xvec(2)-xvec(1);
dy = yvec(2)-yvec(1);

%% AA & bb
AA = sparse([],[],[],ntot,ntot,5*ntot); % alloco lo spazio
tnoto = zeros(ntot,1);

for ii = 2:nx-1
    for jj = 2:ny-1
        kk = ii+nx*(jj-1);
        ke = kk +1;
        kw = kk -1;
        kn = kk +nx;
        ks = kk -nx;
        
        AA(kk,kk) = -2*cond*(1/dx^2+1/dy^2);
        AA(kk,kn) = cond/dy^2;
        AA(kk,ks) = cond/dy^2;
        AA(kk,ke) = cond/dx^2;
        AA(kk,kw) = cond/dx^2;
    end
end

%% coco
% south
jj = 1;
for ii = 1:nx
    kk = ii;
    AA(kk,kk) = 1;
    tnoto(kk) = Ts;
end

% north
jj = ny;
for ii = 1:nx
    kk = ii+nx*(jj-1);
    AA(kk,kk) = 1;
    tnoto(kk) = Tn;
end

% east
ii = nx;
for jj = 1:ny
    kk = ii+nx*(jj-1);
    AA(kk,kk) = 1;
    tnoto(kk) = Te;
end

% west
ii = 1;
for jj = 1:ny
    kk = ii+nx*(jj-1);
    AA(kk,kk) = 1;
    tnoto(kk) = Tw;
end

%% risolvo e surf
TT = AA\tnoto;

[xmat, ymat] = meshgrid(xvec, yvec);
TTmatrix = reshape(TT,nx,ny);

figure(1)
surf(xmat, ymat, TTmatrix', 'facecolor', 'interp')
xlabel('x', 'fontsize', 18)
ylabel('y', 'fontsize', 18)
zlabel('Temperatura', 'fontsize', 18)