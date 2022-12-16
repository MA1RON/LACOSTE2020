clc; clear; close all;
%% dati
LL = 1;
T0 = 500;
T1 = 400;
Teast = 300;
hheast = 10;
cond = .72;

%% griglia spaziale
nx = 51;
ny = 51;
ntot = nx*ny;
xvet = linspace(0,LL,nx);
yvet = linspace(0,LL,ny);
dx = xvet(2)-xvet(1);
dy = yvet(2)-yvet(1);
[xmat,ymat] = meshgrid(xvet,yvet);

%% creo il problema
AA = sparse([],[],[],ntot,ntot,5*ntot);
bb = zeros(ntot,1);

% --- centro ---
for jx = 2:nx-1
    for jy = 2:ny-1
        kk = jx+(jy-1)*nx;
        AA(kk,kk) = -2*(1/dx^2+1/dy^2);
        AA(kk,kk+nx) = 1/dy^2; % north
        AA(kk,kk-nx) = 1/dy^2; % south
        AA(kk,kk+1) = 1/dx^2; % west
        AA(kk,kk-1) = 1/dx^2; % east
    end
end

% --- condizioni al contorno ---
% north (con north-west e north-east) - dirichlet
for jx = 1:nx
    kk = jx+(ny-1)*nx;
    AA(kk,kk) = 1;
    bb(kk) = T0;
end

% south (con south-west e south-east) - dirichlet
for jx = 1:nx
    kk = jx;
    AA(kk,kk) = 1;
    bb(kk) = T0;
end

% west (senza north-west e south-west) - dirichlet
for jy = 2:ny-1
    kk = 1+(jy-1)*nx;
    AA(kk,kk) = 1;
    bb(kk) = T0;
end

% east (senza north-east e south-east) - robin
for jy = 2:ny-1
    kk = jy*nx;
    AA(kk,kk) = -cond*(dx/dy+dy/dx)-hheast*dy;
    AA(kk,kk+nx) = cond*dx/2/dy; % north
    AA(kk,kk-nx) = cond*dx/2/dy; % south
    AA(kk,kk-1) = cond*dy/dx; % west
    bb(kk) = -hheast*Teast*dy;
end

%% --- risolvo condizione iniziale ---
TT = AA\bb;

%% parte il transitorio
tend = 2e5;
dt = tend/100;
time = 0:dt:tend;

ro = 1920;
cp = 835;

Teastave = zeros(length(time),1);

% video settings
loops = length(time)-1;
v = VideoWriter('ex2.avi');
open(v)
M(loops) = struct('cdata',[],'colormap',[]);

inst = 1;
while inst<length(time)
    TTmat = reshape(TT,nx,ny);
    Teastave(inst) = mean(TTmat(nx,:));
    Tp = TT;
    
    II = spdiags(ones(ntot,1),0,ntot,ntot);
    BB = II - AA*dt*cond/ro/cp;
    dd = Tp + dt*bb;
    
    % --- condizioni al contorno ---
    % north (con north-west e north-east) - dirichlet
    for jx = 1:nx
        kk = jx+(ny-1)*nx;
        BB(kk,kk) = 1;
        dd(kk) = T1;
    end

    % south (con south-west e south-east) - dirichlet
    for jx = 1:nx
        kk = jx;
        BB(kk,kk) = 1;
        dd(kk) = T1;
    end

    % west (senza north-west e south-west) - dirichlet
    for jy = 2:ny-1
        kk = 1+(jy-1)*nx;
        BB(kk,kk) = 1;
        dd(kk) = T1;
    end

    % east (senza north-east e south-east) - robin
    for jy = 2:ny-1
        kk = jy*nx;
        BB(kk,kk) = -cond*(dx/dy+dy/dx)-hheast*dy;
        BB(kk,kk+nx) = cond*dx/2/dy; % north
        BB(kk,kk-nx) = cond*dx/2/dy; % south
        BB(kk,kk-1) = cond*dy/dx; % west
        dd(kk) = -hheast*Teast*dy;
    end
    
    % --- risolvo ---
    TT = BB\dd;
    
    % video settings
    surf(xmat,ymat,TTmat','facecolor','interp')
    view([30 15])
    xlabel('x (m)')
    ylabel('y (m)')
    zlabel('Temperature (^oC)')
    hold on
    zlim([300 500])
    
    ax = gca;
    ax.NextPlot = 'replaceChildren';
    drawnow
    M(inst) = getframe;
    
    
    inst = inst +1;
end

% video settings
writeVideo(v,M)
close(v)

%% post production
figure
plot(time(1:inst-1)/3600,Teastave(1:inst-1),'b-','linewidth',2)
hold on
plot([time(1) time(inst-1)]/3600, [Teast Teast],'k-.','linewidth',2)
legend('Media parete est', 'Aria su parete est')
xlabel('Tempo [h]')
ylabel('Temperatura [K]')
ylim([Teast-10 max(Teastave)+10])
grid on
box on
set(gca,'fontsize',18)