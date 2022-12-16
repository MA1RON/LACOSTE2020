clc; clear; close all;
%% dati
% topologia
LL = 3;
DD = 20;
Lqqq = .4;

surf = pi*DD^2/4;
volqqq = Lqqq*surf;

% fisici
ro = 997;
cp = 4182;
RR = 1.5;
II = 50;
GG = .1;
uu = GG/surf/ro;
T0 = 300;

qvol = RR*II^2/volqqq;

%% griglia temporale
t_start = 1;
t_end = 5e6; % guess
dt = (t_end-t_start)/100;
time = (t_start:dt:t_end)';

time_to_plot = round(linspace(t_start,t_end,10)/dt+1);
lines_plotted = 0;

%% griglia spaziale
dx = LL/100;
xx = (0:dx:LL)';
nn = length(xx);

qqq = (xx>LL/2-Lqqq & xx<LL/2+Lqqq).*qvol;

%% risolvo con BE - upwind
TT = ones(nn,1)*T0;

aa = uu*dt/dx;

diag_i = -aa*ones(nn,1);
diag_p = (1+aa)*ones(nn,1);
AA = spdiags([diag_i diag_p],-1:0,nn,nn);
bb = qqq/ro/cp;

AA(1,1) = 1; AA(1,2) = 0;

for inst = 1:length(time)
    Tp = TT;
    BB = Tp+bb;
    BB(1) = T0;
    
    TT = AA\BB;
    if inst == time_to_plot(lines_plotted+1)
        txt_plot = ['t = ',num2str(time(inst)),'s'];
        plot(xx,TT-273.15,'-','Color',[inst/time_to_plot(end) 0 1-inst/time_to_plot(end)],'linewidth',2,'DisplayName',txt_plot)
        if not(inst==length(time_to_plot))
            hold on
            lines_plotted = lines_plotted+1;
        end
        set(gca,'fontsize',16)
        grid on
        box on
        ylabel('Temperatura [^oC]')
        xlabel('Lunghezza [m]')
        legend('location','northwest')
    end
end