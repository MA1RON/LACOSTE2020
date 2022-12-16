clc; clear; close all;
%% proprietà
qqq = -1e5;
cond = 1e4;
ro = 80;
cp = 5;

Lx = 3;
Ly = 1;
Lz = 1;

Tup = 0;
Tdown = 0;
Tnorth = 0;
Tsouth = 0;
qqeast = 0;
hhwest = 1000;
Twest0 = 1000;
T0 = 200;

%% griglia temporale
mm = 100;
time_end = 100;
time = linspace(0,time_end,mm);

dt = time(2)-time(1);

% Twest = Twest0*sin(time*pi*2/time_end);
Twest = linspace(-Twest0,Twest0,mm);

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
TT = ones(ntot,1)*T0; % parto omogeneo
Tp = TT;

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

II = spdiags(ones(ntot,1),0,ntot,ntot);
BB = II-dt*cond/ro/cp*AA;
dd = bb/ro/cp;
TTmat = zeros(nn,nn,nn);

% --- condizioni al contorno ---
% up (con spigoli) - dirichlet
for jx = 1:nn
    for jy = 1:nn
        kk = jx + (jy-1)*nn + (nn-1)*nn^2;
        
        BB(kk,kk) = 1;
        dd(kk) = Tup;
    end
end

% down (con spigoli) - dirichlet
for jx = 1:nn
    for jy = 1:nn
        kk = jx + (jy-1)*nn;
        
        BB(kk,kk) = 1;
        dd(kk) = Tdown;
    end
end

% north (con spigoli west-east) - dirichlet
for jx = 1:nn
    for jz = 2:nn-1
        kk = jx + (nn-1)*nn + (jz-1)*nn^2;
        
        BB(kk,kk) = 1;
        dd(kk) = Tnorth;
    end
end

% south (con spigoli west-east) - dirichlet
for jx = 1:nn
    for jz = 2:nn-1
        kk = jx + (jz-1)*nn^2;
        
        BB(kk,kk) = 1;
        dd(kk) = Tsouth;
    end
end

% east (senza spigoli) - neumann non omogeneo
for jy = 2:nn-1
    for jz = 2:nn-1
        kk = nn + (jy-1)*nn + (jz-1)*nn^2;
        
        BB(kk,kk) = -(dx*dy/dz+dz*dx/dy)-dy*dz/dx*2;
        BB(kk,kk+nn^2) = dx/2*dy/dz; % up
        BB(kk,kk-nn^2) = dx/2*dy/dz; % down
        BB(kk,kk+nn) = dx/2*dz/dy; % north
        BB(kk,kk-nn) = dx/2*dz/dy; % south
        BB(kk,kk-1) = dy*dz/dx*2; % west
        dd(kk) = - qqeast/cond;
    end
end

% west (senza spigoli) - robin con Twest iniziale
for jy = 2:nn-1
    for jz = 2:nn-1
        kk = 1 + (jy-1)*nn + (jz-1)*nn^2;

        BB(kk,kk) = -(dx*dy/dz+dz*dx/dy)-dy*dz/dx*2-hhwest/cond*dy*dz/dx*2;
        BB(kk,kk+nn^2) = dx/2*dy/dz; % up
        BB(kk,kk-nn^2) = dx/2*dy/dz; % down
        BB(kk,kk+nn) = dx/2*dz/dy; % north
        BB(kk,kk-nn) = dx/2*dz/dy; % south
        BB(kk,kk+1) = dy*dz/dx*2; % east
        dd(kk) = -Twest(1)*hhwest/cond*dy*dz/dx*2;
    end
end

%% transitorio
figure
set(gcf,'position',[200,100,2000,700])
subplot(1,2,2)
plot(time,Twest(1:mm),'b-.','linewidth',2)
title('Twest (convezione)')
xlabel('Tempo [s]')
hold on

for instant = 1:mm
    Tp = TT;
    
    dd = Tp+dt*bb;
    
    % --- condizioni al contorno ---
    % up (con spigoli) - dirichlet
    for jx = 1:nn
        for jy = 1:nn
            kk = jx + (jy-1)*nn + (nn-1)*nn^2;
            dd(kk) = Tup;
        end
    end

    % down (con spigoli) - dirichlet
    for jx = 1:nn
        for jy = 1:nn
            kk = jx + (jy-1)*nn;
            dd(kk) = Tdown;
        end
    end

    % north (con spigoli west-east) - dirichlet
    for jx = 1:nn
        for jz = 2:nn-1
            kk = jx + (nn-1)*nn + (jz-1)*nn^2;
            dd(kk) = Tnorth;
        end
    end

    % south (con spigoli west-east) - dirichlet
    for jx = 1:nn
        for jz = 2:nn-1
            kk = jx + (jz-1)*nn^2;
            dd(kk) = Tsouth;
        end
    end

    % east (senza spigoli) - neumann non omogeneo
    for jy = 2:nn-1
        for jz = 2:nn-1
            kk = nn + (jy-1)*nn + (jz-1)*nn^2;
            dd(kk) = - qqeast/cond;
        end
    end

    % west (senza spigoli) - robin con Twest funzione del tempo
    for jy = 2:nn-1
        for jz = 2:nn-1
            kk = 1 + (jy-1)*nn + (jz-1)*nn^2;
            dd(kk) = -Twest(instant)*hhwest/cond*dy*dz/dx*2;
        end
    end
    
    %% risolvo il sistema
    TT = BB\dd;
    
    TTmat = reshape(TT,nn,nn,nn);
    
    %% plotto
    subplot(1,2,1)
    h = slice(xmat,ymat,zmat,TTmat,xvet,yvet,zvet);
    set(h,'EdgeColor','none','FaceAlpha',.07)
    colorbar
    title('Stazionario 3D')
    xlabel('x direction [m]')
    ylabel('y direction [m]')
    zlabel('z direction [m]')
    grid on 
    box on
    set(gca,'fontsize',16)
    drawnow
    
    subplot(1,2,2)
    plot(time(instant),Twest(instant),'ko','markersize',10,'markerfacecolor','y')
end