clc; clear; close all;
%% dati
% geometrici
DD = 4e-3;
LL = 2;

% fisici
ro = 997;
cp = 4187;
kk = .6;

T0x = 20;

% parametri
uu = 5; % m/s
uuli = 0; uuls = 10;

qqq = .1e6; % W/m^3
qqqli = -1e6; qqqls = 1e6;

Tt0 = 10; % °C
Tt0li = 0; Tt0ls = 100;

kk = .6;
kkli = .05; kkls = 1000;

%% dominio
% temporale
dt = 1e-3;

% assiale
dx = LL/200;
xx = (0:dx:LL);
nn = length(xx);

%% figure
% uifigure
f = uifigure('Position',[50 50 1350 800]);

% pannello per slider
pnl = uipanel(f,'Position',[950 50 400 500]);

sTin = uislider(pnl,'Orientation','vertical','Position',[10 10 3 330]);
sTin.Limits = [Tt0li Tt0ls];
sTin.Value = T0x;
lTin = uilabel(pnl,'Position',[5 370 100 20],'Text','T in [°C] ','FontSize',12);
tbTinvalue = uitextarea(pnl,'Position',[5 350 60 20],'Value',num2str(Tt0),'FontSize',14);

sqqq = uislider(pnl,'Orientation','vertical','Position',[85 10 3 330]);
sqqq.Limits = [qqqli qqqls];
sqqq.Value = qqq;
lqqq = uilabel(pnl,'Position',[80 370 100 20],'Text','Q [MW/m^3]','FontSize',12);
tbqqqvalue = uitextarea(pnl,'Position',[80 350 60 20],'Value',num2str(qqq),'FontSize',14);

suuu = uislider(pnl,'Orientation','vertical','Position',[160 10 3 330]);
suuu.Limits = [uuli uuls];
suuu.Value = uu;
luuu = uilabel(pnl,'Position',[155 370 100 20],'Text','U [m/s]','FontSize',12);
tbuuuvalue = uitextarea(pnl,'Position',[155 350 60 20],'Value',num2str(uu),'FontSize',14);

skkk = uislider(pnl,'Orientation','vertical','Position',[235 10 3 330]);
skkk.Limits = [kkli kkls];
skkk.Value = kk;
lkkk = uilabel(pnl,'Position',[230 370 100 20],'Text','k [W/mK]','FontSize',12);
tbkkkvalue = uitextarea(pnl,'Position',[230 350 60 20],'Value',num2str(kk),'FontSize',14);

% grafico
ax = uiaxes(f,'Position',[30 30 900 770],...
            'XLim',[0 LL],'YLim',[0 100],...
            'FontSize',18,'Box','on',...
            'XGrid','on','YGrid','on');
ax.YLabel.String = 'T [^oC]';
ax.XLabel.String = 'x [m]';
ax.Title.String = 'CONDUCTION PLUS ADVECTION';

% bottoni per uscire
bg = uibuttongroup(f,'Position',[950 650 200 70]);
cstay = uiradiobutton(bg,'Position',[90 5 100 35],'Text','Stai','Value',true);
cleave = uiradiobutton(bg,'Position',[10 5 100 35],'Text','Chiudi','Value',false);

% label proprietà
bgprop = uibuttongroup(f,'Position',[950 710 200 70]);

lro = uilabel(bgprop,'Position',[10 35 40 20],'Text','ro = ','FontSize',14);
tbro = uitextarea(bgprop,'Position',[52 35 50 20],'Value',int2str(ro),'FontSize',14);
lroum = uilabel(bgprop,'Position',[104 35 50 20],'Text','kg/m^3','FontSize',14);
lcp = uilabel(bgprop,'Position',[10 5 40 20],'Text','cp = ','FontSize',14);
tbcp = uitextarea(bgprop,'Position',[52 5 50 20],'Value',int2str(cp),'FontSize',14);
lcpum = uilabel(bgprop,'Position',[104 5 50 20],'Text','J/kgK','FontSize',14);

% bottoni per oscillare
bgoscil = uibuttongroup(f,'Position',[950 550 200 100]);
csino = uilistbox(bgoscil,'Position',[30 5 130 90],...
                  'Items',{'-','T','Q','u'},...
                  'Multiselect','on');

%% BE - upwind
TT = ones(nn,1)*T0x;
Tp = TT;
dT = (Tt0ls-Tt0li)/100;
dq = (qqqls-qqqli)/100;
du = (uuls-uuli)/100;
isgoingup = true;

while not(cleave.Value)
    %% oscillo
    if any(char(csino.Value) == 'T') && isgoingup
        sTin.Value = sTin.Value+dT;
        if sTin.Value >= Tt0ls - .1*(Tt0ls-Tt0li)
            isgoingup = false;
        end
    elseif any(char(csino.Value) == 'T') && not(isgoingup)
        sTin.Value = sTin.Value-dT;
        if sTin.Value <= Tt0li + .1*(Tt0ls-Tt0li)
            isgoingup = true;
        end
    end
    if any(char(csino.Value) == 'Q') && isgoingup
        sqqq.Value = sqqq.Value+dq;
        if sqqq.Value >= qqqls - .1*(qqqls-qqqli)
            isgoingup = false;
        end
    elseif any(char(csino.Value) == 'Q') && not(isgoingup)
        sqqq.Value = sqqq.Value-dq;
        if sqqq.Value <= qqqli + .1*(qqqls-qqqli)
            isgoingup = true;
        end
    end
    if any(char(csino.Value) == 'u') && isgoingup
        suuu.Value = suuu.Value+du;
        if suuu.Value >= uuls - .1*(uuls-uuli)
            isgoingup = false;
        end
    elseif any(char(csino.Value) == 'u') && not(isgoingup)
        suuu.Value = suuu.Value-du;
        if suuu.Value <= uuli + .1*(uuls-uuli)
            isgoingup = true;
        end
    end
    
    %% cambio valori
    if str2double(tbTinvalue.Value) ~= Tt0 && (str2double(tbTinvalue.Value) >= Tt0li && str2double(tbTinvalue.Value) <= Tt0ls)
        Tt0 = str2double(tbTinvalue.Value);
        sTin.Value = Tt0;
    elseif sTin.Value ~= Tt0
        Tt0 = sTin.Value;
        tbTinvalue.Value = num2str(Tt0);
    end
    if str2double(tbqqqvalue.Value) ~= qqq && (str2double(tbqqqvalue.Value) >= qqqli && str2double(tbqqqvalue.Value) <= qqqls)
        qqq = str2double(tbqqqvalue.Value);
        sqqq.Value = qqq;
    elseif sqqq.Value ~= qqq
        qqq = sqqq.Value;
        tbqqqvalue.Value = num2str(qqq);
    end
    if str2double(tbuuuvalue.Value) ~= uu && (str2double(tbuuuvalue.Value) >= uuli && str2double(tbuuuvalue.Value) <= uuls)
        uu = str2double(tbuuuvalue.Value);
        suuu.Value = uu;
    elseif suuu.Value ~= uu
        uu = suuu.Value;
        tbuuuvalue.Value = num2str(uu);
    end
    if str2double(tbkkkvalue.Value) ~= kk && (str2double(tbkkkvalue.Value) >= kkli && str2double(tbkkkvalue.Value) <= kkls)
        kk = str2double(tbkkkvalue.Value);
        skkk.Value = kk;
    elseif skkk.Value ~= kk
        kk = skkk.Value;
        tbkkkvalue.Value = num2str(kk);
    end
    
    ro = str2double(tbro.Value(1)); % kg/m^3
    cp = str2double(tbcp.Value(1)); % J/kgK
    %% risolvo il sistema
    
    diag_i = (-uu*dt/dx-kk*dt/dx^2/ro/cp)*ones(nn,1);
    diag_p = (1+uu*dt/dx+2*kk*dt/dx^2/ro/cp)*ones(nn,1);
    diag_s = (-kk*dt/dx^2/ro/cp)*ones(nn,1);
    AA = spdiags([diag_i diag_p diag_s],-1:1,nn,nn);
    bb = Tp + qqq/ro/cp;
    
    % coco
    AA(1,1) = 1;
    AA(1,2) = 0;
    bb(1) = Tt0;
    
    AA(end,end) = 1+uu*dt/dx;
    AA(end,end-1) = -uu*dt/dx;
    bb(end) = Tp(end) + qqq/ro(end)/cp(end);
    
    TT = AA\bb; Tp = TT;
    
    plot(ax,xx,TT,'r-','linewidth',2);
    drawnow
end

close(f)