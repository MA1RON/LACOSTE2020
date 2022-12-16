clc; clear; close all;
%% dati
LL = .05;

sigma =  5.670374419e-8;
Tsun = @(t) (400)*abs(sin(pi*t/24/60));

Tair = @(t) 10*abs(sin(pi*(t-1)/24/60));
hhext = 50;
Tint = 20;
hhint = 25;

ro = 2400;
cp = 1000;
kk = 1.13;
alpha = kk/ro/cp;

%% griglia t
tend = 3*24*60;
dt = tend/500;
time = (0:dt:tend)';
mm = length(time);

%% griglia x
dx = LL/100;
xx = (0:dx:LL)';
nn = length(xx);

%% BE
optimset.Display = 'off';
TT = 12.5*ones(nn,1); Tp = TT;

for inst = 2:mm
    fun = @(T) funbuilder(T,Tp,dx,dt,sigma,kk,alpha,Tsun(time(inst)),Tair(time(inst)),Tint,hhint,hhext);
    TT = fsolve(@(T) fun(T),Tp,optimset);
    Tp = TT;
    
    subplot(1,2,1)
    plot(xx,TT,'b-','linewidth',2)
    set(gca,'fontsize',18)
    grid on
    box on
    xlabel('x [m]')
    ylabel('T [^oC]')
    ylim([0 25])
    xlim([0 LL])
    
    subplot(1,2,2)
    plot(time(inst)/60,Tair(time(inst)),'ko','markerfacecolor','b','markersize',10)
    set(gca,'fontsize',18)
    grid on
    box on
    xlabel('tempo [h]')
    ylabel('T [^oC]')
    ylim([0 25])
    xlim([0 tend/60])
    drawnow
end

%% funzioni
function fun = funbuilder(TT,Tp,dx,dt,sigma,kk,alpha,Tsun,Tair,Tint,hhint,hhext)
nn = length(TT);

fun(1) = TT(2)-TT(1)+sigma/kk*dx*(Tsun^4-TT(1)^4)+hhext/kk*dx*(Tair-TT(1));
fun(2:nn-1) = TT(2:end-1)-Tp(2:end-1)-alpha*dt*(TT(3:end)-2*TT(2:end-1)+TT(1:end-2))/dx^2;
fun(nn) = TT(end-1)-TT(end)+hhint/kk*dx*(Tint-TT(end));
end