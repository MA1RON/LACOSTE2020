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

%% creo il problema ( AA & bb )
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
T0 = AA\bb;

%% griglia temporale
tend = 24;
dt = 1/60;
time = 0:dt:tend;

%% parte il transitorio
ro = 7.800; % dato presumibilmente sbagliato da 7800
cp = 450;
Text = @(t) 5 + 15*sin(2*pi*t/24); % t in ore
TT = T0;
TTcenter = zeros(length(time),1);

inst = 1;
while inst<length(time)
    TTcenter(inst) = TT(nx1*(ny1-1)+round(nx1/2));
    
    Tp = TT;
    
    II = spdiags(ones(ntot,1),0,ntot,ntot);
    BB = II-dt*cond/ro/cp*AA;
    dd = Tp+dt*bb;
    
    % --- condizioni al contorno ---
    % north (con north-west and north-east)
    for jx = 1:nx2
        kk = nx1*(ny1-1)+nx2*(ny2-ny1)+jx;
        BB(kk,kk) = 1;
        dd(kk) = Text(time(inst));
    end

    % south esterno
    for jx = 2:nx1-1
        kk = jx;
        BB(kk,kk) = -2*(1/dx^2+1/dy^2);
        BB(kk,kk+nx1) = 2/dy^2; % north
        BB(kk,kk-1) = 1/dx^2; % west
        BB(kk,kk+1) = 1/dx^2; % east
        % dd(kk) = 0 neumann omogeneo
    end

    % south interno (con nodo centrale and east-center)
    for jx = nx1:nx2
        kk = nx1*(ny1-1)+jx;
        BB(kk,kk) = 1;
        dd(kk) = Tint;
    end

    % west (con south-west)
    for jy = 1:ny2-1
        if jy<=ny1
            kk = (jy-1)*nx1+1;
        else
            kk = (ny1-1)*nx1+(jy-ny1)*nx2+1;
        end
        BB(kk,kk) = 1;
        dd(kk) = Text(time(inst));
    end

    % east esterno
    for jy = ny1+1:ny2-1
        kk = nx1*(ny1-1)+(jy-ny1+1)*nx2;
        BB(kk,kk) = -2*(1/dx^2+1/dy^2);
        BB(kk,kk+nx2) = 1/dy^2; % north
        BB(kk,kk-nx2) = 1/dy^2; % south
        BB(kk,kk-1) = 2/dx^2; % west
        % dd(kk) = 0 neumann omogeneo per simmetria
    end

    % east interno (con south-center)
    for jy = 1:ny1-1
        kk = nx1+(jy-1)*nx1;
        BB(kk,kk) = 1;
        dd(kk) = Tint;
    end
    
    
    TT = BB\dd;
    
    %{
    %% reshape & plot ogni ora
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
    if mod(inst,60) == 0
        surf(xmat,ymat,TTmat','facecolor','interp')
        pause(2)
        close
    end
    %}
    
    inst = inst +1;
end

%% post production - transitorio in un punto
figure(1)
plot(time(1:inst-1),TTcenter(1:inst-1),'b-','linewidth',2)
hold on
plot(time(1:inst-1),Text(time(1:inst-1)),'r--','linewidth',2)
grid on
box on
xlabel('Tempo [h]')
ylabel('Temperatura [^oC]')
set(gca,'fontsize',18)
legend('Temperatura punto centrale', 'Temperatura esterna', 'location', 'northeastoutside')

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
            TTmat(jx,jy) = T0(kk);
        end
    end
end

%% post production - stazionario
figure(2)
surf(xmat,ymat,TTmat','facecolor','interp')
xlabel('Direzione x [m]')
ylabel('Direzione y [m]')
zlabel('Temperatura [^oC]')
set(gca,'fontsize',18)
title('Stazionario 2D')
close(2)