clc; clear; close all;
%% dati
Lx = 27e-3/2; % simmetrico
Ly = 12e-3;

condmain = 5;
condchip = 50;
qqqchip = 1e7;
hh = 500;
Tinf = 20;

Tlim = 85;

%% griglia spaziale
nx2 = 51;
ny2 = 51;
nx1 = round(2/3*nx2);
ny1 = round(3/4*ny2);
ntot = nx2*ny2;
xvet = linspace(0,Lx,nx2);
yvet = linspace(0,Ly,ny2);
dx = xvet(2)-xvet(1);
dy = yvet(2)-yvet(1);

[xmat, ymat] = meshgrid(xvet,yvet);

%% creo il problema
AA = sparse([],[],[],ntot,ntot,5*ntot);
bb = zeros(ntot,1);

% --- centro ---
for jx = 2:nx2-1
    for jy = 2:ny2-1
        kk = jx+(jy-1)*nx2;
        
        AA(kk,kk) = -2*(1/dx^2+1/dy^2);
        AA(kk,kk+nx2) = 1/dy^2; % north
        AA(kk,kk-nx2) = 1/dy^2; % south
        AA(kk,kk-1) = 1/dx^2; % west
        AA(kk,kk+1) = 1/dx^2; % east
        
        % zona chip (nx1 & ny1 esclusi)
        if jx>nx1 && jy>ny1
            bb(kk) = -qqqchip/condchip;
        end
    end
end

% --- bordo di attacco ---
% orizzonatale
for jx = nx1+1:nx2-1
    kk = jx+(ny1-1)*nx2;
    
    AA(kk,kk) = -(condchip+condmain)*dx/dy-(condmain+condchip)*dy/dx;
    AA(kk,kk+nx2) = dx/dy*condchip; % north
    AA(kk,kk-nx2) = dx/dy*condmain; % south
    AA(kk,kk-1) = dy/dx/2*(condchip+condmain); % west
    AA(kk,kk+1) = dy/dx/2*(condchip+condmain); % east
end

% verticale
for jy = ny1+1:ny2-1
    kk = nx1+(jy-1)*nx2;
    
    AA(kk,kk) = -(condchip+condmain)*dy/dx-(condmain+condchip)*dx/dy;
    AA(kk,kk+nx2) = dy/dx/2*(condchip+condmain); % north
    AA(kk,kk-nx2) = dy/dx/2*(condchip+condmain); % south
    AA(kk,kk-1) = dx/dy*condmain; % west
    AA(kk,kk+1) = dx/dy*condchip; % east
end

% angolo
for jx = nx1:nx1
    kk = jx+(ny1-1)*nx2;
    
    AA(kk,kk) = -(condmain+condchip)/2*(dx/dy+dy/dx)-condmain*dx/dy-condmain*dy/dx;
    AA(kk,kk+nx2) = (condmain+condchip)/2*dx/dy; % north
    AA(kk,kk-nx2) = condmain*dx/dy; % south
    AA(kk,kk-1) = condmain*dy/dx; % west
    AA(kk,kk+1) = (condmain+condchip)/2*dy/dx; % east
end

% --- condizioni al contorno ---
%{
% north dirichlet
for jx = 1:nx2-1
    kk = (ny2-1)*nx2+jx;
    
    AA(kk,kk) = 1;
    bb(kk) = Tinf;
end
%}

% north-west (robin north neumann west)
for jx = 1:1
    kk = jx+(ny2-1)*nx2;
    
    AA(kk,kk) = -condmain/2*(dx/dy+dy/dx)-dx/2*hh;
    AA(kk,kk-nx2) = condmain*dx/2/dy; % south
    AA(kk,kk+1) = condmain*dy/2/dx; % east
    bb(kk) = -dx/2*hh*Tinf;
end


% north (senza north-west e north-east)
for jx = 2:nx2-1
    kk = jx+(ny2-1)*nx2;
    if jx>nx1
        cond = condchip;
    else
        cond = condmain;
    end
    
    AA(kk,kk) = -dx/dy*cond-dy/dx*cond-hh*dx;
    AA(kk,kk-nx2) = dx/dy*cond; % south
    AA(kk,kk-1) = dy/dx/2*cond; % west
    AA(kk,kk+1) = dy/dx/2*cond; % east
    bb(kk) = -dx*hh*Tinf;
    
    %{
    if jx == nx1+1 % cambio west
        AA(kk,kk-1) = dy/dx/2*condmain; % west <- !
    elseif jx == nx1 % cambio east
        AA(kk,kk+1) = dy/dx/2*condchip; % east <- !
    end
    %}
end

% north-east (robin north neumann east per simmetria)
for jx = nx2:nx2
    kk = jx+(ny2-1)*nx2;
    
    AA(kk,kk) = -condchip/2*(dx/dy+dy/dx)-dx/2*hh;
    AA(kk,kk-nx2) = condchip*dx/2/dy; % south
    AA(kk,kk-1) = condchip*dy/2/dx; % west
    bb(kk) = -dx/2*hh*Tinf;
end

% east (senza north-east e south-east)
for jy = 2:ny2-1
    kk = jy*nx2;
    if jy>ny1
        cond = condchip;
    else
        cond = condmain;
    end
    
    AA(kk,kk) = -cond*(dx/dy+dy/dx);
    AA(kk,kk+nx2) = cond*dx/2/dy; % north
    AA(kk,kk-nx2) = cond*dx/2/dy; % south
    AA(kk,kk-1) = cond*dy/dx; % west
    
    %{
    if jy == ny1+1
        AA(kk,kk-nx2) = condmain*dx/2/dy; % south
    elseif jy == ny1
        AA(kk,kk+nx2) = condchip*dx/2/dy; % north
    end
    %}
end

%{
% east dirichlet
for jy = 2:ny2-1
    kk = jy*nx2;
    
    AA(kk,kk) = 1;
    bb(kk) = Tinf;
end
%}

% south-east (nuemann south neumann east)
for jy = 1:1
    kk = nx2;
    
    AA(kk,kk) = -(dx/dy+dy/dx)/2;
    AA(kk,kk+nx2) = dx/2/dy; % north
    AA(kk,kk-1) = dy/2/dx; % west
end

% south (senza south-east e south-west)
for jx = 2:nx2-1
    kk = jx;
    
    AA(kk,kk) = -(dx/dy+dy/dx);
    AA(kk,kk+nx2) = dx/dy; % north
    AA(kk,kk-1) = dy/2/dx; % west
    AA(kk,kk+1) = dy/2/dx; % east
end

% south-west (neumann south e neumann west)
for jx = 1:1
    kk = jx;
    
    AA(kk,kk) = -dx/dy/2-dy/dx/2;
    AA(kk,kk+nx2) = dx/2/dy; % north
    AA(kk,kk+1) = dy/2/dx; % west
end

% west (senza north-west e south-west)
for jy = 2:ny2-1
    kk = 1+(jy-1)*nx2;
    
    AA(kk,kk) = -(dx/dy+dy/dx);
    AA(kk,kk+nx2) = dx/2/dy; % north
    AA(kk,kk-nx2) = dx/2/dy; % south
    AA(kk,kk+1) = dy/dx; % east
end

% --- risolvo ---
TT = AA\bb;

%% verifico temperatura massima
Tmax = max(TT);
if Tmax>Tlim
    fprintf('La temperatura massima raggiunta supera il limite.\n')
else
    fprintf('La temperatura massima raggiunta non supera il limite.\n')
end
fprintf('T massima: %.3f °C. T limite: %.3f °C.\n', Tmax, Tlim)

%% post production
TTmat = reshape(TT,nx2,ny2);

surf(xmat,ymat,TTmat','facecolor','interp')
xlabel('x [m]')
ylabel('y [m]')
zlabel('T [^oC]')
set(gca,'fontsize',18)
title('Stazionario 2D')
colorbar