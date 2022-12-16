clc; clear; close all;
%% dati
% geometrici
DD = 4e-3;
LL = 2;

% fisici
ro = 997;
cp = 4187;

T0x = 20;

% parametri
uu = 5; % m/s
uuli = 0; uuls = 50;

qqq = 0; % W/m^3
qqqli = -1e6; qqqls = 1e6;

Tt0 = 20; % °C
Tt0li = 0; Tt0ls = 100;

kk = 0;
kkli = 0; kkls = 100;

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
tbTinvalueli = uitextarea(pnl,'Position',[5 400 60 20],'Value',num2str(Tt0li),'FontSize',14);
tbTinvaluels = uitextarea(pnl,'Position',[5 430 60 20],'Value',num2str(Tt0ls),'FontSize',14);
lTinvl = uilabel(pnl,'Position',[5 450 100 20],'Text','T in lim [°C] ','FontSize',12);

sqqq = uislider(pnl,'Orientation','vertical','Position',[85 10 3 330]);
sqqq.Limits = [qqqli qqqls];
sqqq.Value = qqq;
lqqq = uilabel(pnl,'Position',[80 370 100 20],'Text','Q [W/m^3]','FontSize',12);
tbqqqvalue = uitextarea(pnl,'Position',[80 350 60 20],'Value',num2str(qqq),'FontSize',14);
tbqqqvalueli = uitextarea(pnl,'Position',[80 400 60 20],'Value',num2str(qqqli),'FontSize',14);
tbqqqvaluels = uitextarea(pnl,'Position',[80 430 60 20],'Value',num2str(qqqls),'FontSize',14);
lqqqvl = uilabel(pnl,'Position',[75 450 100 20],'Text','Q lim [W/m^3] ','FontSize',12);

suuu = uislider(pnl,'Orientation','vertical','Position',[160 10 3 330]);
suuu.Limits = [uuli uuls];
suuu.Value = uu;
luuu = uilabel(pnl,'Position',[155 370 100 20],'Text','U [m/s]','FontSize',12);
tbuuuvalue = uitextarea(pnl,'Position',[155 350 60 20],'Value',num2str(uu),'FontSize',14);
tbuuuvalueli = uitextarea(pnl,'Position',[155 400 60 20],'Value',num2str(uuli),'FontSize',14);
tbuuuvaluels = uitextarea(pnl,'Position',[155 430 60 20],'Value',num2str(uuls),'FontSize',14);
luuuvl = uilabel(pnl,'Position',[155 450 100 20],'Text','U lim [m/s] ','FontSize',12);

skkk = uislider(pnl,'Orientation','vertical','Position',[235 10 3 330]);
skkk.Limits = [kkli kkls];
skkk.Value = kk;
lkkk = uilabel(pnl,'Position',[230 370 100 20],'Text','k [W/mK]','FontSize',12);
tbkkkvalue = uitextarea(pnl,'Position',[230 350 60 20],'Value',num2str(kk),'FontSize',14);
tbkkkvalueli = uitextarea(pnl,'Position',[230 400 60 20],'Value',num2str(kkli),'FontSize',14);
tbkkkvaluels = uitextarea(pnl,'Position',[230 430 60 20],'Value',num2str(kkls),'FontSize',14);
lkkkvl = uilabel(pnl,'Position',[230 450 100 20],'Text','k lim [W/mK] ','FontSize',12);

% grafico
ax = uiaxes(f,'Position',[30 30 900 770],...
            'XLim',[0 LL],'YLim',[0 100],...
            'FontSize',18,'Box','on',...
            'XGrid','on','YGrid','on');
ax.YLabel.String = 'T [^oC]';
ax.XLabel.String = 'x [m]';
ax.Title.String = 'CONDUCTION PLUS ADVECTION';

% tasto per mettere punto
bgpoint = uibutton(f,'state','Position',[950 660 200 40],'Text','Inserisci punto');

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
              
% tasto per uscire
bgexit = uibutton(f,'state','Position',[1150 710 200 70],'Text','Esci');


%% BE - upwind
TT = ones(nn,1)*T0x;
Tp = TT;

isgoingup = true;
n_point = 0;

max_single_point = 10;
x_points = zeros(max_single_point,1);

while not(bgexit.Value)
    %% punto fisso
    % lo creo
    if bgpoint.Value && max_single_point > n_point
        bgpoint.Value = false;
        
        n_point = n_point + 1;
        x_points(n_point) = 0;
    end
    % lo mostro
    for p = 1:n_point
        x_points(p) = x_points(p) + uu*dt/dx;
        j_points(p) = ceil(x_points(p)/LL);
        
        if j_points(p) >= nn
            n_point = n_point - 1;
            j_points = j_points(2:end);
            x_points = x_points(2:end);
            if n_point == 0
                hold(ax,'off')
            end
            break
        end
        
        if p == 1
            hold(ax,'off')
        end
        plot(ax,xx(j_points(p)),TT(j_points(p)),'ko','markersize',10,'markerfacecolor','y')
        hold(ax,'on')
    end
    
    %% oscillo
    if any(char(csino.Value) ~= '-')
        if any(char(csino.Value) == 'T') && isgoingup
            dT = (Tt0ls-Tt0li)/100;
            sTin.Value = sTin.Value+dT;
            if sTin.Value >= Tt0ls - .1*(Tt0ls-Tt0li)
                isgoingup = false;
            end
        elseif any(char(csino.Value) == 'T') && not(isgoingup)
            dT = (Tt0ls-Tt0li)/100;
            sTin.Value = sTin.Value-dT;
            if sTin.Value <= Tt0li + .1*(Tt0ls-Tt0li)
                isgoingup = true;
            end
        end
        if any(char(csino.Value) == 'Q') && isgoingup
            dq = (qqqls-qqqli)/100;
            sqqq.Value = sqqq.Value+dq;
            if sqqq.Value >= qqqls - .1*(qqqls-qqqli)
                isgoingup = false;
            end
        elseif any(char(csino.Value) == 'Q') && not(isgoingup)
            dq = (qqqls-qqqli)/100;
            sqqq.Value = sqqq.Value-dq;
            if sqqq.Value <= qqqli + .1*(qqqls-qqqli)
                isgoingup = true;
            end
        end
        if any(char(csino.Value) == 'u') && isgoingup
            du = (uuls-uuli)/100;
            suuu.Value = suuu.Value+du;
            if suuu.Value >= uuls - .1*(uuls-uuli)
                isgoingup = false;
            end
        elseif any(char(csino.Value) == 'u') && not(isgoingup)
            du = (uuls-uuli)/100;
            suuu.Value = suuu.Value-du;
            if suuu.Value <= uuli + .1*(uuls-uuli)
                isgoingup = true;
            end
        end
    end
    
    %% controllo limiti
    if (char(csino.Value) == '-')
        if str2double(tbTinvalueli.Value(1)) ~= Tt0li && Tt0 >= str2double(tbTinvalueli.Value(1))
            Tt0li = str2double(tbTinvalueli.Value(1));
            sTin.Limits = [Tt0li Tt0ls];
            tbTinvalueli.Value = num2str(Tt0li);
        end
        if str2double(tbTinvaluels.Value(1)) ~= Tt0ls && Tt0 <= str2double(tbTinvaluels.Value(1))
            Tt0ls = str2double(tbTinvaluels.Value(1));
            sTin.Limits = [Tt0li Tt0ls];
            tbTinvaluels.Value = num2str(Tt0ls);
        end
        if str2double(tbqqqvalueli.Value(1)) ~= qqqli && qqq >= str2double(tbqqqvalueli.Value(1))
            qqqli = str2double(tbqqqvalueli.Value(1));
            sqqq.Limits = [qqqli qqqls];
            tbqqqvalueli.Value = num2str(qqqli);
        end
        if str2double(tbqqqvaluels.Value(1)) ~= qqqls && qqq <= str2double(tbqqqvaluels.Value(1))
            qqqls = str2double(tbqqqvaluels.Value(1));
            sqqq.Limits = [qqqli qqqls];
            tbqqqvaluels.Value = num2str(qqqls);
        end
        if str2double(tbuuuvalueli.Value(1)) ~= uuli && uu >= str2double(tbuuuvalueli.Value(1)) && str2double(tbuuuvalueli.Value(1)) >= 0
            uuli = str2double(tbuuuvalueli.Value(1));
            suuu.Limits = [uuli uuls];
            tbuuuvalueli.Value = num2str(uuli);
        end
        if str2double(tbuuuvaluels.Value(1)) ~= uuls && uu <= str2double(tbuuuvaluels.Value(1))
            uuls = str2double(tbuuuvaluels.Value(1));
            suuu.Limits = [uuli uuls];
            tbuuuvaluels.Value = num2str(uuls);
        end
        if str2double(tbkkkvalueli.Value(1)) ~= kkli && kk >= str2double(tbkkkvalueli.Value(1))
            kkli = str2double(tbkkkvalueli.Value(1));
            skkk.Limits = [kkli kkls];
            tbkkkvalueli.Value = num2str(kkli);
        end
        if str2double(tbkkkvaluels.Value(1)) ~= kkls && kk <= str2double(tbkkkvaluels.Value(1))
            kkls = str2double(tbkkkvaluels.Value(1));
            skkk.Limits = [kkli kkls];
            tbkkkvaluels.Value = num2str(kkls);
        end
    end
    
    %% cambio valori
    if str2double(tbTinvalue.Value(1)) ~= Tt0 && (str2double(tbTinvalue.Value(1)) >= Tt0li && str2double(tbTinvalue.Value(1)) <= Tt0ls)
        Tt0 = str2double(tbTinvalue.Value(1));
        sTin.Value = Tt0;
        tbTinvalue.Value = num2str(Tt0);
    elseif sTin.Value ~= Tt0
        Tt0 = sTin.Value;
        tbTinvalue.Value = num2str(Tt0);
    end
    if str2double(tbqqqvalue.Value(1)) ~= qqq && (str2double(tbqqqvalue.Value(1)) >= qqqli && str2double(tbqqqvalue.Value(1)) <= qqqls)
        qqq = str2double(tbqqqvalue.Value(1));
        sqqq.Value = qqq;
        tbqqqvalue.Value = num2str(qqq);
    elseif sqqq.Value ~= qqq
        qqq = sqqq.Value;
        tbqqqvalue.Value = num2str(qqq);
    end
    if str2double(tbuuuvalue.Value(1)) ~= uu && (str2double(tbuuuvalue.Value(1)) >= uuli && str2double(tbuuuvalue.Value(1)) <= uuls)
        uu = str2double(tbuuuvalue.Value(1));
        suuu.Value = uu;
        tbuuuvalue.Value = num2str(uu);
    elseif suuu.Value ~= uu
        uu = suuu.Value;
        tbuuuvalue.Value = num2str(uu);
    end
    if str2double(tbkkkvalue.Value(1)) ~= kk && (str2double(tbkkkvalue.Value(1)) >= kkli && str2double(tbkkkvalue.Value(1)) <= kkls)
        kk = str2double(tbkkkvalue.Value(1));
        skkk.Value = kk;
        tbkkkvalue.Value = num2str(kk);
    elseif skkk.Value ~= kk
        kk = skkk.Value;
        tbkkkvalue.Value = num2str(kk);
    end
    
    if str2double(tbro.Value(1)) ~= ro
        ro = str2double(tbro.Value(1)); % kg/m^3
        tbro.Value = num2str(ro);
    end
    if str2double(tbcp.Value(1)) ~= cp
        cp = str2double(tbcp.Value(1)); % J/kgK
        tbcp.Value = num2str(cp);
    end
    
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