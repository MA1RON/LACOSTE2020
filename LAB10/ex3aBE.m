clc; clear; close all;
%% dati
% geometrici
DD = 5e-3;
LL = 1;

surf = pi*DD^2/4;
vol = surf*LL;
Awet = pi*DD*LL;

% fisici
ro = 140;
cp = 6000;
uu = 1e-2;
T0 = 5;
Tleft = 10;
THe = 5;
hh = 100;
time_to_plot = [5 10 25 50 100];
line_plotted = 0;

%% griglia temporale
tend = time_to_plot(end);
dt = tend/100;
time = (0:dt:tend)';

%% griglia spaziale
dx = LL/100;
xx = (0:dx:LL)';
nn = length(xx);

%% risolvo con BE - upwind
TT = ones(nn,1)*T0;

aa = uu*dt/dx;
diag_i = -aa*ones(nn,1);
diag_p = (1+aa+hh*Awet/vol*dt/ro/cp)*ones(nn,1);
AA = spdiags([diag_i diag_p],-1:0,nn,nn);
bb = hh*Awet/vol/ro/cp*THe*dt*ones(nn,1);

AA(1,1) = 1; AA(1,2) = 0;

jj = 2;
while jj<=length(time)
    Tp = TT;
    BB = Tp+bb;
    BB(1) = Tleft;
    
    TT = AA\BB;
    
    if time(jj) >= time_to_plot(line_plotted+1)
        plot(xx*100,TT,'linewidth',2,'Color',[jj/length(time) 0 1-jj/length(time)])
        if jj == length(time)
            legend('\Delta=5s','\Delta=10s','\Delta=25s','\Delta=50s','\Delta=100s')
        else
            hold on
        end
        grid on
        box on
        set(gca,'fontsize',18)
        xlabel('Lunghezza [cm]')
        ylabel('Temperatura [K]')
        title('Transitorio advection')
        
        line_plotted = line_plotted+1;
    end
    
    jj = jj+1;
end