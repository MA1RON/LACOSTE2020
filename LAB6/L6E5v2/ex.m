clc; clear; close all;
%% dati
Tacqua = 70+273.15;
hh = 10e3;
LL = 2e-2;
kk = 400;
Tplasma = 2000;
Tlim = 550;

sigma = 1.380649e-23;

%% griglia
dx = 1e-4;
xx = (0:dx:LL)';
nn = length(xx);

%% soluzione con newton
TTguess = 1500*ones(nn,1);
tol = 1e-3;
jjmax = 1e3;

fun = @(T) funbuilder(T,dx,sigma,kk,Tplasma,hh,Tacqua);
jac = @(T) jacbuilder(T,dx,sigma,kk,Tplasma,hh,Tacqua);

tic
[TT,res,err,lastjj] = newton(fun,jac,TTguess,tol,jjmax);
timenewana = toc;

fprintf('NEWTON ANALITICO:\tjj = %i\terr = %.3e\tres = %.3e.\n',lastjj,err,res)

plot(xx,TT,'b-.','linewidth',2,'DisplayName','Newton - analitiche')


%% funzioni utili
function fun = funbuilder(TT,dx,sigma,kk,Tplasma,hh,Tacqua)
    fun = zeros(size(TT));
    
    fun(1) = TT(2)-TT(1) + dx*sigma/kk*(Tplasma^4-TT(1)^4);
    fun(2:end-1) = TT(3:end)-2*TT(2:end-1)+TT(1:end-2);
    fun(end) = TT(end-1)-TT(end) - dx*hh/kk*(TT(end)-Tacqua);
end

function jac = jacbuilder(TT,dx,sigma,kk,Tplasma,hh,Tacqua)
    diag_i = zeros(size(TT));
    diag_p = zeros(size(TT));
    diag_s = zeros(size(TT));
    
    diag_i(1:end-2) = 1;
    diag_i(end-1) = 1;
    diag_i(end) = 0;
    
    diag_p(1) = 1;
    diag_p(2:end-1) = -1 + dx*sigma/kk*-4*TT(1)^3;
    diag_p(end) = -1-dx*hh/kk;
    
    diag_s(1) = 0;
    diag_s(2) = 1;
    diag_s(3:end) = 1;
    
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