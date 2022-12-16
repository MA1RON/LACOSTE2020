clc; clear; close all;
%% dati
% topologia
LL = 100;
DD = 11e-3;
tt = 4e-2;

surf = pi*DD^2/4;
vol = LL*surf;
Awet = (7.5e-2)^2;

% fisici
kksolido = 1;
kkfluido = .63;
ro = 994;
cp = 4186;
ni = 719e-6;
QQ = 3*60/1000;
uu = QQ/surf;
T0 = 20;
Tleft = 35;

% hh_int
Rein = ro*uu*DD/ni;
Prin = cp*ni/kksolido;
Nuin = .3+(.62*Rein^.5*Prin^(1/3)*(1+(.4/Prin)^(2/3))^-.25)*(1+(Rein/282000)^.625)^.8;
hhin = Nuin*kkfluido/DD;

Rsupint = .1;
Rcond = tt/kksolido;

Raou = 0; % guess no wind
Nuou = (.671+.37*Raou^(1/6))^2;
hhou = Nuou*kkfluido/LL;

UU = 1/(1/hhin+Rsupint+Rcond+1/hhou);

%% griglia temporale
t_end = .1; % guess
dt = t_end/100;
time = (0:dt:t_end)';

time_to_plot = round(linspace(0,t_end,10)/dt+1);
lines_plotted = 0;

%% griglia spaziale
dx = LL/100;
xx = (0:dx:LL)';
nn = length(xx);

%% risolvo con BE - upwind
TT = ones(nn,1)*T0;
finalT = zeros(length(time),1);

aa = uu*dt/dx;
gamma = UU*Awet/vol/ro/cp*dt;

diag_i = -aa*ones(nn,1);
diag_p = (1+aa+gamma)*ones(nn,1);
AA = spdiags([diag_i diag_p],-1:0,nn,nn);
bb = gamma*T0*ones(nn,1);

AA(1,1) = 1; AA(1,2) = 0;

tol = 1e-5;
inst = 1;
err = 1;
while inst<length(time) && err>tol
    Tp = TT;
    BB = Tp+bb;
    BB(1) = Tleft;
    
    TT = AA\BB;
    if inst == time_to_plot(lines_plotted+1)
        txt_plot = ['t = ',num2str(time(inst)),'s'];
        plot(xx,TT,'-','Color',[inst/time_to_plot(end) 0 1-inst/time_to_plot(end)],'linewidth',2,'DisplayName',txt_plot)
        if not(inst==length(time_to_plot))
            hold on
            lines_plotted = lines_plotted+1;
        end
        set(gca,'fontsize',16)
        grid on
        box on
        ylabel('Temperatura [^oC]')
        xlabel('Lunghezza [m]')
        legend('location','northeastoutside')
        title('Transitorio nel tubo')
    end
    
    finalT(inst) = TT(end);
    err = norm(TT-Tp)/norm(TT-T0);
    inst = inst + 1;
end
if err<=tol
    fprintf('Tolleranza raggiunta al tempo t=%.3f\n',time(inst))
else
    fprintf('Tolleranza non raggiunta')
end

%% plotto T finale
figure
plot(time(1:inst-1),finalT(1:inst-1),'b-','linewidth',2)
set(gca,'fontsize',16)
grid on
box on
ylabel('Temperatura [^oC]')
xlabel('Tempo [s]')
legend('T_{finale}','location','southeast')
title('Transitorio temperatura finale del tubo')