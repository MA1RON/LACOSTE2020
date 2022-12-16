clc; clear; close all;
%% dati
% geometrici
LL = 100;
HH = 5.96;
Din = 64e-3;
Dou = 70e-3;

surfpiastra = LL*HH;
surftubo = pi/4*(Dou^2-Din^2);
voltubo = surftubo*LL;

% fisici
% percNaNO3 = .6; % dato inutile (non usato)
% percKNO3 = .4; % dato inutile (non usato)
qq = 800;
% qlin = 4768; % dato inutile (qq*surfpiastra/LL)
GG = 1.2;
T0 = 290;
cp = 1500;
ro = @(T) (4.444e-07*T.^2-1.389e-3*T+2.519)*1e3;

qqq = qq*surfpiastra/voltubo;
uu = @(T) GG./ro(T)/surftubo;

%% plotto dati per avere un idea di cosa succederà
figure
set(gca,'fontsize',18)
TT_guess = linspace(T0,5e3);
yyaxis left
plot(TT_guess,ro(TT_guess),'b-','linewidth',2)
ylabel('Densità [kg/m^3]')
hold on
yyaxis right
plot(TT_guess,uu(TT_guess),'r-','linewidth',2)
legend('Densità','Velocità', 'location', 'north')
xlabel('Temperatura [^oC]')
ylabel('Velocità [m/s]')
title('Dati funzione della Temperatura')
grid on
box on

%% griglia spaziale
dx = LL/100;
xx = (0:dx:LL)';
nn = length(xx);

%% griglia temporale
tend = 100; % guess
dt = tend/100;
time = (0:dt:tend)';
mm = length(time);

%% risolvo con BE - upwind
figure
howmanylines = 7;
time_to_plot = round(linspace(0,tend/dt,howmanylines)+1);
lines_plotted = 0;

TT = ones(nn,1) * T0;

aa = @(T) uu(T)*dt/dx;

inst = 1;
while inst<=mm
    Tp = TT;
    
    % frozen coefficients
    diag_i = -aa(Tp);
    diag_p = (1+aa(Tp));
    AA = spdiags([diag_i diag_p],-1:0,nn,nn);
    bb = qqq./ro(Tp)/cp*dt;
    
    AA(1,1) = 1; AA(1,2) = 0;
    BB = Tp+bb; BB(1) = T0;
    
    TT = AA\BB;
    
    if inst == time_to_plot(lines_plotted+1)
        text = ['t = ',int2str(time(inst)),'s'];
        plot(xx,TT,'-','linewidth',2,'Color',[inst/mm 0 1-inst/mm], 'DisplayName',text)
        hold on
        legend('location', 'northwest')
        xlabel('Lunghezza [m]')
        ylabel('Temperatura [^oC]')
        title('Transitorio pure advection')
        set(gca,'fontsize',18)
        grid on
        box on
        
        lines_plotted = lines_plotted+1;
    end
    
    inst = inst+1;
end
