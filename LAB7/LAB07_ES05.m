%% LaCoSTe 2020/2021 
%% LAB07_ES05 - Soluzioni, Matteo Nicoli
close all
clear all
clc
 
% La soluzione analitica del problema proposto (dT(t)/dt=sin(t)+T) è
% T(t)=1/2*exp(t)-1/2*(cos(t)+sin(t))
 
%Soluzione analitica valutata all'istante di tempo finale (per calcolo
%errori)
an=(1/2*exp(1)-1/2*(cos(1)+sin(1)));
 
dt=logspace(-1,-3,3)';
t0=0; t1=1;
ff=@(t)(sin(t));
color=['r';'g';'b';'c'];
T0=0;

err1=zeros(size(dt));
err2=zeros(size(dt));
err3=zeros(size(dt));
h1=zeros(size(dt));
h2=zeros(size(dt));
h3=zeros(size(dt));
 
for ii=1:length(dt)
    tt=t0:dt(ii):t1;
    
    %Risoluzione del problema con i 3 metodi richiesti, implementati nelle
    %rispettive funzioni esterne.
    [T1,time1]=ES05_BE(t0,t1,T0,ff,dt(ii));
    [T2,time2]=ES05_FE(t0,t1,T0,ff,dt(ii));
    [T3,time3]=ES05_CN(t0,t1,T0,ff,dt(ii));
    
    err1(ii)=abs(T1(end)-an)/an;
    err2(ii)=abs(T2(end)-an)/an;
    err3(ii)=abs(T3(end)-an)/an;
    
    figure(1)
    hold on;
    h1(ii)=plot(tt,T1,color(ii,:),'linewidth',1.5);
    figure(2)
    hold on;
    h2(ii)=plot(tt,T2,color(ii,:),'linewidth',1.5);
    figure(3)
    hold on;
    h3(ii)=plot(tt,T3,color(ii,:),'linewidth',1.5);
end

% Risoluzione del problema con le funzioni Matlab "ode23" e "ode45".
fun = @(t,T)(sin(t)+T);
tspan = [0 1];
[tout23,Tout23] = ode23(fun,tspan,T0);
[tout45,Tout45] = ode45(fun,tspan,T0);
 
%% POST PROCESSING

leg = [{['\Deltat=',num2str(dt(1))]},...
       {['\Deltat=',num2str(dt(2))]},...
       {['\Deltat=',num2str(dt(3))]},{'Analitica'}];

tan=t0:0.05:t1;
Tan=(1/2*exp(tan)-1/2*(cos(tan)+sin(tan)));
figure(1)
set(gca,'fontsize',16)
title('Backward Euler');
grid minor
box on
xlabel('t')
ylabel('T')
h1(ii+1)=plot(tan,Tan,'ko','linewidth',1.5);
legend(h1,leg,'location','nw');
print -depsc ES05FE
 
figure(2)
set(gca,'fontsize',16)
title('Forward Euler');
grid minor
box on
xlabel('t')
ylabel('T')
h2(ii+1)=plot(tan,Tan,'ko','linewidth',1.5);
legend(h2,leg,'location','nw');
print -depsc ES05BE
 
figure(3)
set(gca,'fontsize',16)
title('Crank-Nicolson');
grid minor
box on
xlabel('t')
ylabel('T')
h3(ii+1)=plot(tan,Tan,'ko','linewidth',1.5);
legend(h3,leg,'location','nw');
print -depsc ES05CN
 
figure(4)
set(gca,'fontsize',16)
loglog(dt,err1,'r-*',dt,err2,'b-o',dt,err3,'k-s','linewidth',2,'markersize',12)
grid minor
xlabel('\Deltat (s)')
ylabel('Errore relativo')
title('Convergenza')
legend('Backward Euler','Forward Euler','Crank-Nicolson','location','SE')
print -depsc -f1 ES05converg

figure(5)
set(gca,'fontsize',16)
title('Confronto con funzioni MATLAB');
hold on
grid minor
box on
xlabel('t')
ylabel('T')
h4(1)=plot(tt,T1,'b-s','linewidth',3,'markersize',3,'markerfacecolor','b'); % FE
h4(2)=plot(tt,T2,'r-d','linewidth',2.5,'markersize',3,'markerfacecolor','r'); % BE
h4(3)=plot(tt,T3,'g-o','linewidth',2,'markersize',3,'markerfacecolor','g'); % CN
h4(4)=plot(tout23,Tout23,'m-.','linewidth',1.5); % ode23
h4(5)=plot(tout45,Tout45,'c--','linewidth',1); % ode45
h4(6)=plot(tan,Tan,'k-o','linewidth',1); % analytical
legend(h4,'FE','BE','CN','ode23','ode45','Analitica','location','nw');
print -depsc ES05ode
