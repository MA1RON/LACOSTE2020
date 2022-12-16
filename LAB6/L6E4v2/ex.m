clc; clear; close all;
%% dati
LL = 1;
DD = .2;

Awet = pi*DD*LL;
vol = pi*DD^2/4*LL;

TTinf = 300+273.15;
hh = 50;

tol = 1e-2;

qqq = @(x) 8e3*cos(pi*x/LL);
kk = @(T) 100./(11.8+.0238*(T)) + 8.775e-11*(T).^3;

%% griglia
dx = 1e-3;
xx = (-LL/2:dx:0)'; % dominio simmetrico
nn = length(xx);

TTguess = TTinf*ones(nn,1);

%% plotting
figure
xlabel('Lunghezza [m]')
ylabel('Temperatura [°C]')
set(gca,'fontsize',18)
grid on
box on
hold on
legend('location','northwest')

%% NEWTON DERIVATE ANALITICHE
jjmax = 1e3;

dkkdT = @(T)(-0.0238*100./(11.8+0.0238*(T)).^2+3*8.775e-11*(T).^2);
d2kkdT2 = @(T)(2*0.0238^2*100./(11.8+0.0238*(T)).^3+6*8.775e-11*(T));

fun = @(T) funbuilder(T,kk,dkkdT,qqq(xx),TTinf,hh,dx,Awet,vol);
jac = @(T) jacbuilder(T,kk,dkkdT,d2kkdT2,qqq(xx),TTinf,hh,dx,Awet,vol);

tic
[TT,res,err,lastjj] = newton(fun,jac,TTguess,tol,jjmax);
timenewana = toc;

fprintf('NEWTON ANALITICO:\tjj = %i\terr = %.3e\tres = %.3e.\n',lastjj,err,res)

plot(xx,TT,'b-.','linewidth',2,'DisplayName','Newton - analitiche')

%% NEWTON DERIVATE NUMERICHE
jjmax = 1e3;

dT = tol;
dkkdT = @(T) (kk(T+dT)-kk(T-dT))/2/dT;
d2kkdT2 = @(T) (dkkdT(T+dT)-dkkdT(T-dT))/2/dT;

fun = @(T) funbuilder(T,kk,dkkdT,qqq(xx),TTinf,hh,dx,Awet,vol);
jac = @(T) jacbuilder(T,kk,dkkdT,d2kkdT2,qqq(xx),TTinf,hh,dx,Awet,vol);

tic
[TT,res,err,lastjj] = newton(fun,jac,TTguess,tol,jjmax);
timenewnum = toc;

fprintf('NEWTON NUMERICO:\tjj = %i\terr = %.3e\tres = %.3e.\n',lastjj,err,res)

plot(xx,TT,'r-.','linewidth',2,'DisplayName','Newton - numeriche')

%% MATLAB fsolve (numerico)dT = tol;
dkkdT = @(T)(-0.0238*100./(11.8+0.0238*(T)).^2+3*8.775e-11*(T).^2);

fun = @(T) funbuilder(T,kk,dkkdT,qqq(xx),TTinf,hh,dx,Awet,vol);

OPTIONS = optimoptions('fsolve','OptimalityTolerance',tol,'Display','off');
tic
[TT,res,~,OUTPUT,~] = fsolve(fun,TTguess,OPTIONS);
timematfso = toc;

fprintf('MATLAB fsolve:\t\tjj = %i\tres = %.3e.\n',OUTPUT.iterations,norm(res))

plot(xx,TT,'k-.','linewidth',2,'DisplayName','fsolve')

%% iterazioni di punto fisso
jjmax = 1e3;

tic
[TT,err,lastjj] = fixedpoint(kk,qqq(xx),dx,hh,TTinf,TTguess,tol,jjmax,Awet,vol,TTinf);
timefixpoi = toc;

fprintf('PUNTO FISSO:\t\tjj = %i\terr = %.3e\n',lastjj,err)

plot(xx,TT,'g-.','linewidth',2,'DisplayName','Punto fisso')

%% confronto
figure
bar(categorical({'Newton analitiche','Newton numeriche','fsolve','Fixed point'}),[timenewana,timenewnum,timematfso,timefixpoi])
ylabel('Tempo [s]')
set(gca,'yscale','log')
set(gca,'fontsize',18)
grid on
box on

%% funzioni
function fun = funbuilder(TT,kk,dkkdT,qvet,TTinf,hh,dx,Awet,vol)
    fun = zeros(size(TT));
    
    fun(1) = TT(2)-TT(1) - hh/kk(TT(1))*dx*(TT(1)-TTinf);
    fun(2:end-1) = dkkdT(TT(2:end-1)).*((TT(3:end)-TT(1:end-2))/(2*dx)).^2 + kk(TT(2:end-1)).*(TT(1:end-2)-2*TT(2:end-1)+TT(3:end))/dx^2 + qvet(2:end-1) + hh*(TT(2:end-1)-TTinf)*Awet/vol;
    fun(end) = TT(end) - TT(end-1);
end

function jac = jacbuilder(TT,kk,dkkdT,d2kkdT2,qvet,TTinf,hh,dx,Awet,vol)
    diag_s = ones(size(TT));
    diag_p = ones(size(TT));
    diag_i = ones(size(TT));
    
    diag_s(1) = 0;
    diag_s(2) = 1;
    diag_s(3:end) = dkkdT(TT(2:end-1)).*(TT(3:end)-TT(1:end-2))/2/dx^2 + kk(TT(2:end-1))/dx^2;
    
    diag_p(1) = -1-hh/kk(TT(1))*dx;
    diag_p(2:end-1) = d2kkdT2(TT(2:end-1)).*(TT(3:end)-TT(1:end-2))/4/dx^2 + dkkdT(TT(2:end-1)).*(TT(1:end-2)-2*TT(2:end-1)+TT(3:end))/dx^2 + kk(TT(2:end-1))*-2/dx^2 + hh*Awet/vol;
    diag_p(end) = 1;
    
    diag_i(1:end-2) = -dkkdT(TT(2:end-1)).*(TT(3:end)-TT(1:end-2))/2/dx^2 + kk(TT(2:end-1))/dx^2;
    diag_i(end-1) = -1;
    diag_i(end) = 0;
    
    jac = spdiags([diag_i diag_p diag_s],-1:1,length(TT),length(TT));
end

function [TT,res,err,jj] = newton(fun,jac,TTguess,tol,jjmax)
    err = 1; res = 1;
    TT = TTguess;
    Tp = TTguess;
    
    jj = 1;
    while err > tol || res > tol
        
        TT = TT - jac(TT)\fun(TT);
        
        err = norm(TT-Tp)/norm(TT-TTguess);
        res = norm(fun(TT));
        
        if jj >= jjmax
            fprintf('Newton non converge.\n')
            return
        end
        
        Tp = TT;
        jj = jj + 1;
    end
end

function [TT,err,jj] = fixedpoint(kk,qvet,dx,hh,Tinf,TTguess,tol,jjmax,Awet,vol,TTinf)
    TT = TTguess;
    Tp = TTguess;
    
    err = 1;
    
    jj = 1;
    while err > tol
        % primi k ma sbagliati: occhio!
        % kkminus = (kk([TT(2:end);TT(end)])-kk([TT(1);TT(1:end-1)]))/4 + kk(TT(1:end));
        % kkplus = -(kk([TT(2:end);TT(end)])-kk([TT(1);TT(1:end-1)]))/4 + kk(TT(1:end));
        
        kkminus = (4*kk(TT)+kk([TT(1);TT(1:end-1)])-kk([TT(2:end);TT(end)]))/4;
        kkplus  = (4*kk(TT)-kk([TT(1);TT(1:end-1)])+kk([TT(2:end);TT(end)]))/4;
        % kkave = (kkminus+kkplus)/2;
        
        diag_s = [0;kkplus(1:end-1)];
        diag_p = -kkminus-kkplus + hh*Awet/vol*dx^2;
        diag_i = [kkminus(2:end);0];
        
        AA = spdiags([diag_i diag_p diag_s],-1:1,length(TTguess),length(TTguess));
        bb = - (qvet + hh*TTinf*Awet/vol)*dx^2;
        
        % --- condizioni al contorno ---
        % robin
        AA(1,1) = -1-dx*hh/kk(TT(1));
        AA(1,2) = 1;
        bb(1) = dx*hh/kk(TT(1))*-Tinf;
        
        % neumann omogeneo
        AA(end,end-1) = 1;
        AA(end,end) = -1;
        bb(end) = 0; 
        
        TT = AA\bb;
        
        err = norm(TT-Tp)/norm(TT-TTguess);
        
        Tp = TT;
        
        if jj >= jjmax
            fprintf('Punto fisso non converge.\n')
            return
        end
        
        jj = jj + 1;
    end
end