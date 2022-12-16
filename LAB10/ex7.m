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
qq = @(t) -8.36e-7*t.^2+2.24e-2*t+664;
% qlin = 4768; % dato inutile (qq*surfpiastra/LL)
GG = 1.2;
T0 = 290;
cp = 1500;
ro = @(T) (4.444e-07*T.^2-1.389e-3*T+2.519)*1e3;

qqq = @(t) qq(t)*surfpiastra/voltubo;
uu = @(T) GG./ro(T)/surftubo;

%% griglia temporale
tend = (18-8)*3600;
dt = tend/100;
time = (0:dt:tend)';
mm = length(time);

%% plotto dati per avere un idea di cosa succederà
% uu and ro
figure
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
set(gca,'fontsize',18)

% qq
figure
plot(time/3600+8,qq(time),'b-','linewidth',2)
legend('Irradianza', 'location', 'south')
xlabel('Tempo [h]')
ylabel('Irradianza [W/m^2]')
title('Dato funzione del Tempo')
grid on
box on
set(gca,'fontsize',18)

%% griglia spaziale
dx = LL/100;
xx = (0:dx:LL)';
nn = length(xx);

%% risolvo con BE - upwind
figure
howmanylines = 13;
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
    bb = qqq(time(inst))./ro(Tp)/cp*dt;
    
    AA(1,1) = 1; AA(1,2) = 0;
    BB = Tp+bb; BB(1) = T0;
    
    TT = AA\BB;
    
    if inst == time_to_plot(lines_plotted+1)
        hour = floor(time(inst)/3600+8);
        minute = (ceil(time(inst)/3600)-time(inst)/3600)*60;
        if length(int2str(minute)) == 1
            minute = [int2str(minute),'0'];
        else
            minute = int2str(minute);
        end
        text = ['t = ',int2str(hour),':',minute];
        plot(xx,TT,'-.','linewidth',2,'Color',[inst/mm 0 1-inst/mm], 'DisplayName',text)
        hold on
        legend('location', 'northwestoutside')
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
