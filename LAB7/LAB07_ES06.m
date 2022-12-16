%% LaCoSTe 2020/2021 
%% LAB07_ES06 - Soluzioni, Matteo Nicoli
clear all
close all
clc
 
% Traiettoria esatta.
figure(1);
set(gca,'fontsize',16)
tt=0:0.01:2*pi;
x_ex=4*cos(tt)+2;
y_ex=3*sin(tt)+3;
plot(x_ex,y_ex,'linewidth',2);
title('Traiettoria analitica');
xlabel('x [-]')
ylabel('y [-]')
axis equal
box on;
grid minor;
print -depsc -f1 ES06Fig1
 
%Definisco diversi numeri di "time step" per la suddivisione dell'intervallo di tempo.
nt=logspace(1,4,4);
color=['b','g','m','k','r'];
figure('units','normalized','outerposition',[0.25 0 0.5 0.8])
figure(2);
set(gca,'fontsize',16)
for nn=1:length(nt)
    % Time step
    dt=(2*pi)/nt(nn);
    % Preallocazione
    xy=zeros(2,nt(nn)+1);
    % Condizione iniziale
    xy(:,1)=[2;6];
    % Matrice dei coefficienti
    A=[1 (4/3)*dt;(-3/4)*dt 1];
    % Parte costante del secondo termine
    b=[4;-3/2]*dt;
    % Soluzione con Backward Euler
    for ii=2:nt(nn)+1
        xy(:,ii)=A\(b+xy(:,ii-1));
    end
    % Plotting
    h1(nn)=plot(xy(1,:),xy(2,:),color(nn),'linewidth',2);
    hold on
end
xlabel('x [-]')
ylabel('y [-]')
title('Traiettoria numerica')
box on;
grid minor;
axis equal
% legend(ll,'location','northoutside','orientation','horizontal')
legend(h1,'10 steps','100 steps','1000 steps','10000 steps','location','southoutside')
print -depsc -f2 ES06Fig2

%% solution with ode
y0 = [2;6];
tspan = [0 2*pi];
[~,yout] = ode45(@ES06fun,tspan,y0);
figure('units','normalized','outerposition',[0 0 0.5 0.8])
figure(3)
set(gca,'fontsize',16)
plot(x_ex,y_ex,'linewidth',2);
hold on
plot(xy(1,:),xy(2,:),'r-.','linewidth',2)
plot(yout(:,1),yout(:,2),'k--','linewidth',2)
xlabel('x [-]')
ylabel('y [-]')
axis equal
title('Confronto')
legend('Analitica','Backward Euler (10000 steps)','ode45','location','southoutside')
box on;
grid minor;
print -depsc -f3 ES06Fig3

