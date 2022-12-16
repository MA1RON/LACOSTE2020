%% LaCoSTe 2020/2021 
%% LAB07_ES07 - Soluzioni, Matteo Nicoli
clear all
close all
clc

% Parametri del modello Lotka-Volterra
alpha=1;
beta=0.01;
gamma=0.02;
delta=1;
% Condizione iniziale
x0=20;
y0=20;
t_fin=20;
% Diversi numeri di "time steps"
nt=[1000 10000];
for nn=1:length(nt)
    tic
    % Time step
    dt=(t_fin)/nt(nn);
    tt=linspace(0,t_fin,nt(nn)+1);
    xx = zeros(size(tt));
    yy = zeros(size(tt));
    % Condizione iniziale
    xx(1)=x0;
    yy(1)=y0;
    % Soluzione con metodo Forward Euler
    for ii=2:nt(nn)+1
        xx(ii)=xx(ii-1)+dt*(alpha*xx(ii-1)-beta*yy(ii-1)*xx(ii-1));
        yy(ii)=yy(ii-1)+dt*(gamma*xx(ii-1)*yy(ii-1)-delta*yy(ii-1));
    end
    toc
    % Plotting
    figure('units','normalized','outerposition',[0 0 0.7 0.7])
    figure(nn)
    hold on
    grid minor
    box on
    set(gca,'fontsize',16)
    plot(tt,xx,'r','linewidth',2)
    plot(tt,yy,'b','linewidth',2)
    legend('Prede','Predatori','location','NorthWest')
    if nn==1
        fff=text(5,300,'Instabilità \rightarrow La soluzione è divergente!!!');
        set(fff,'fontsize',16)
    end
    title(sprintf('Simulazione con %.0f time steps',nt(nn)))
    xlabel('Tempo')
    ylim([0 500])
	print('-depsc',['-f' num2str(nn)],['ES07Fig' num2str(nn)]);
end

% Soluzione con "ode23"
tic
tspan = [0 t_fin];
z0 = [x0; y0];
[tout,yout] = ode23(@(tt,yy)ES07fun(tt,yy,alpha,beta,gamma,delta),tspan,z0);
toc
figure('units','normalized','outerposition',[0 0 0.7 0.7])
figure(3)
hold on
grid minor
box on
set(gca,'fontsize',16)
plot(tt,xx,'r-s','linewidth',2)
plot(tt,yy,'b-s','linewidth',2)
plot(tout,yout(:,1),'m--*','linewidth',2)
plot(tout,yout(:,2),'c--*','linewidth',2)
xlabel('Tempo')
legend('Prede (FE)','Predatori (FE)','Prede (ode23)','Predatori (ode23)','location','Northoutside','orientation','horizontal')
print -depsc -f3 ES07Fig3

% Soluzione con metodo Backward Euler (a coefficienti congelati)
tic
nt = 1e4;
dt=(t_fin)/nt;
tt=linspace(0,t_fin,nt+1);
zz = zeros(length(tt),2);
zz(1,:)=[x0 y0];
for ii=2:nt+1
    bb = zz(ii-1,:)'; % bb deve essere un vettore colonna
    AA = [1-alpha*dt           beta*dt*zz(ii-1,1);
          -gamma*dt*zz(ii-1,2) 1+delta*dt       ];
    zz(ii,:) = (AA\bb)'; % zz deve essere un vettore riga
end
toc
% Plotting
figure('units','normalized','outerposition',[0 0 0.7 0.7])
figure(4)
hold on
grid minor
box on
set(gca,'fontsize',16)
plot(tt,zz(:,1),'r-','linewidth',2)
plot(tt,zz(:,2),'b-','linewidth',2)
plot(tout,yout(:,1),'m--','linewidth',2)
plot(tout,yout(:,2),'c--','linewidth',2)
xlabel('Tempo')
legend('Prede (BE+FC)','Predatori (BE+FC)','Prede (ode23)','Predatori (ode23)','location','Northoutside','orientation','horizontal')
print -depsc -f4 ES07Fig4