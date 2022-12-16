clc; clear; close all;
%% dati
% topologia
Lx = 3;
Ly = 1;

% fisici
Tn = 10; % dirichlet
Texts = 20; % robin
hhs = 10;
Tw = 30; % dirichlet
qqe = 40; % neumann
cond = 1;
qqq = 1e4;

%% griglia
nx = 20;
ny = 30;
ntot = nx*ny;
xvec = (linspace(0,Lx,nx));
yvec = (linspace(0,Ly,ny));
dx = xvec(2)-xvec(1);
dy = yvec(2)-yvec(1);

%% risolvo
% inizializzo
AA = sparse([],[],[],nx,ny,5*nx*ny);
bb = zeros(ntot,1);

% centro
for ii = 2:nx-1
    for jj = 2:ny-1
        kk = (jj-1)*nx+ii;
        kn = kk-1;
        ks = kk+1;
        kw = kk-nx;
        ke = kk+nx;
        AA(kk,kk) = -2*cond*(1/dx^2+1/dy^2);
        AA(kk,kn) = cond/dx^2;
        AA(kk,ks) = cond/dx^2;
        AA(kk,kw) = cond/dy^2;
        AA(kk,ke) = cond/dy^2;
        bb(kk) = -qqq;
    end
end

%% boundaries
% north
ii = 1;
for jj = 1:ny
    kk = (jj-1)*nx+ii;
    AA(kk,kk) = 1;
    bb(kk) = Tn;
end

% south
ii = nx;
for jj = 1:ny
    kk = (jj-1)*nx+ii;
    if jj == 1 % southwest -> dirichlet
        AA(kk,kk) = 1;
        bb(kk) = Tw;
    elseif jj == ny % southeast -> neumann
        kn = kk-1;
        kw = kk-nx;
        AA(kk,kk) = -(1/dx^2+1/dy^2);
        AA(kk,kn) = 1/dx;
        AA(kk,kw) = 1/dy;
        bb(kk) = -(1/dx+1/dy)*qqe/cond/dx;
    else % -> robin
        kn = kk-1;
        % ks doesn't exist
        kw = kk-nx;
        ke = kk+nx;
        AA(kk,kk) = -2*(1/dx^2+1/dy^2+hhs/cond/dy);
        AA(kk,kn) = 1/dx^2;
        % AA(kk,ks) = 0;
        AA(kk,kw) = 2/dy^2;
        AA(kk,ke) = 1/dx^2;
        bb(kk) = -2*hhs*Texts/cond/dy;
    end
end

% west
jj = 1;
for ii = 1:nx
    kk = (jj-1)*nx+ii;
    AA(kk,kk) = 1;
    bb(kk) = Tw;
end

% east
jj = ny;
for ii = 2:nx-1
    kk = (jj-1)*nx+ii;
    kn = kk-1;
    ks = kk+1;
    kw = kk-nx;
    % ke doesn't exist
    AA(kk,kk) = -2*(1/dx^2+1/dy^2);
    AA(kk,kn) = 1/dx^2;
    AA(kk,ks) = 1/dx^2;
    AA(kk,kw) = 2/dy^2;
    % AA(kk,ke) = 0;
    bb(kk) = -2*qqe/cond/dy;
end

%% risolvo e surf
TT = AA\bb;

[xmat, ymat] = meshgrid(xvec, yvec);
TTmatrix = reshape(TT,nx,ny);

figure(1)
surf(xmat, ymat, TTmatrix', 'facecolor', 'interp')
xlabel('x', 'fontsize', 18)
ylabel('y', 'fontsize', 18)
zlabel('Temperatura', 'fontsize', 18)