clc; clear; close all;
%% dati
% geometrici
DD = 4e-3;
LL = 2;

% fisici
ro = 130;
cp = 5500;
kk = 2000;

T0x = 20;

% parametri
uu = 5; % mm/s
uuli = 1; uuls = 10;

qqq = .1; % MW/m^3
qqqli = -1; qqqls = 1;

Tt0 = 10; % °C
Tt0li = 0; Tt0ls = 100;

%% dominio
% temporale
tend = 500;
dt = tend/50;
time = (0:dt:tend)';
mm = length(time);

% assiale
dx = LL/200;
xx = (0:dx:LL);
nn = length(xx);

%% figure
% uifigure
f = uifigure('Position',[200 200 1200 800]);

% pannello per slider
pnl = uipanel(f,'Position',[950 50 200 500]);

sTin = uislider(pnl,'Orientation','vertical','Position',[10 10 3 430]);
sTin.Limits = [Tt0li Tt0ls];
sTin.Value = T0x;
lTin = uilabel(pnl,'Position',[5 450 100 20],'Text','T in [°C] ','FontSize',12);

sqqq = uislider(pnl,'Orientation','vertical','Position',[70 10 3 430]);
sqqq.Limits = [qqqli qqqls];
sqqq.Value = qqq;
lqqq = uilabel(pnl,'Position',[65 450 100 20],'Text','Q [MW/m^3]','FontSize',12);

suuu = uislider(pnl,'Orientation','vertical','Position',[130 10 3 430]);
suuu.Limits = [uuli uuls];
suuu.Value = uu;
luuu = uilabel(pnl,'Position',[135 450 100 20],'Text','U [mm/s]','FontSize',12);

% grafico
ax = uiaxes(f,'Position',[30 30 900 770],...
            'XLim',[0 LL],'YLim',[0 100],...
            'FontSize',18,'Box','on',...
            'XGrid','on','YGrid','on');
ax.YLabel.String = 'T [^oC]';
ax.XLabel.String = 'x [m]';
ax.Title.String = 'ADVECTION';

% bottoni per uscire
bg = uibuttongroup(f,'Position',[950 600 200 70]);
cstay = uiradiobutton(bg,'Position',[90 5 100 35],'Text','Stai','Value',true);
cleave = uiradiobutton(bg,'Position',[10 5 100 35],'Text','Chiudi','Value',false);

% label proprietà
bgprop = uibuttongroup(f,'Position',[950 710 200 70]);
lro = uilabel(bgprop,'Position',[10 35 200 35],'Text',['ro = ' int2str(ro) ' kg/m^3'],'FontSize',14);
lcp = uilabel(bgprop,'Position',[10 5 200 35],'Text',['cp = ' int2str(cp) ' J/kgK'],'FontSize',14);

%% BE - upwind
TT = ones(nn,1)*T0x;
Tp = TT;

while not(cleave.Value)
    Tt0 = sTin.Value;
    qqq = sqqq.Value*1e6; % MW/m^3
    uu = suuu.Value*1e-3; % mm/s
    
    diag_i = (-uu*dt/dx)*ones(nn,1);
    diag_p = (1+uu*dt/dx)*ones(nn,1);
    AA = spdiags([diag_i diag_p],-1:0,nn,nn);
    bb = Tp + qqq/ro/cp;
    
    % coco
    AA(1,1) = 1;
    AA(1,2) = 0;
    bb(1) = Tt0;
    
    TT = AA\bb; Tp = TT;
    
    plot(ax,xx,TT,'r-','linewidth',2);
    drawnow
end

close(f)