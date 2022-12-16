clc; clear; close all;
%% dati
% geometrici
Lx2=0.1;
Lx1=0.2;
Ly2=0.1;
Ly1=0.05;

% fisici
Text = 5;
Tint = 20;

cond = 10;

%% griglia spaziale
nx2=21;
nx1=11;
ny2=21;
ny1=11;
ntot=nx1*(ny2-ny1)+nx2*ny1;

xvec=linspace(0,Lx1,nx2);
yvec=linspace(0,Ly2,ny2);
dx=xvec(2)-xvec(1);
dy=yvec(2)-yvec(1);
[xmat,ymat]=meshgrid(xvec,yvec);

%% creo il dominio
TT_dominio = zeros(size(xmat));
TT_dominio(1:nx1,1:ny1-1)=Text;
TT_dominio(nx1+1:nx2,1:ny1-1)=nan;
TT_dominio(1:nx2,ny1:ny2)=Text;
% surf(xmat,ymat,TT_dominio') % mostro il dominio

%% creo il sistema
AA=sparse([],[],[],ntot,ntot,5*ntot);
bb=zeros(ntot,1);

% --- center ---
% center south block
for jx=2:nx1-1
    for jy=2:ny1
        kk=jx+(jy-1)*nx1;
        AA(kk,kk) = -2*cond*(1/dx^2+1/dy^2);
        if jy < ny1
            AA(kk,kk+nx1) = cond/dy^2;
        else % jy == ny1
            AA(kk,kk+nx2) = cond/dy^2;
        end
        AA(kk,kk-nx1) = cond/dy^2;
        AA(kk,kk+1) = cond/dx^2;
        AA(kk,kk-1) = cond/dx^2;
    end 
end 

% center black north
for jx=2:nx2-1
   for jy=ny1+1:ny2-1
       kk = jx + (jy-ny1-1)*nx2 + nx1*ny1 + (nx2-nx1);
       AA(kk,kk) = -2*cond*(1/dx^2+1/dy^2);
       AA(kk,kk+nx2) = cond/dy^2;
       AA(kk,kk-nx2) = cond/dy^2;
       AA(kk,kk+1) = cond/dx^2;
       AA(kk,kk-1) = cond/dx^2;
   end
end 
           
% --- condizioni al contorno ---
% north
jy=ny2;
for jx=1:nx2
    kk=(ny1-1)*nx1+(jy-ny1)*nx2+jx;
    AA(kk,kk)=1;
    bb(kk)=Text;
end 

% south esterno
for jx=2:nx1-1
    kk=jx;
    
    AA(kk,kk) = -2*cond*(1/dy^2+1/dx^2);
    AA(kk,kk+nx1) = 2*cond/dy^2;
    AA(kk,kk-1) = cond/dx^2;
    AA(kk,kk+1) = cond/dx^2;
end 

% south interno
for jx=nx1:nx2
    kk = nx1*ny1 - nx1 + jx;
    AA(kk,kk)=1;
    bb(kk)=Tint;
end 

% west
jx=1;
for jy=1:ny2
    if jy <= ny1
        kk=(jy-1)*nx1+jx;
    else
        kk=(jy-ny1)*nx2+(ny1-1)*nx1+jx;
    end
    AA(kk,kk)=1;
    bb(kk)=Text;
end 

% east esterno (di simmmetria)
jx=nx2;
for jy=ny1+1:ny2-1
    kk = ny1*nx1 + (nx2-nx1) + jx*(jy-ny1);  
    
    AA(kk,kk) = -2*cond*(1/dy^2+1/dx^2);
    AA(kk,kk+nx2) = cond/dy^2;
    AA(kk,kk-nx2) = cond/dy^2;
    AA(kk,kk-1) = 2*cond/dx^2;
end 

% east interno
jx = nx1;
for jy=1:ny1-1
    kk=jx+(jy-1)*nx1;
    AA(kk,kk)=1;
    bb(kk)=Tint;
end 

% --- risolvo ---
TT = AA\bb;

% --- reshape ---
TTmat = zeros(nx1,ny1);
for jx = 1:nx2
    for jy = 1:ny2
        if jx>nx1 && jy<ny1
            TTmat(jx,jy) = nan;
        elseif jy <= ny1
            kk = jx+(jy-1)*nx1;
            TTmat(jx,jy) = TT(kk);
        else
            kk = jx+(jy-1-ny1)*nx2+nx1*ny1 + (nx2-nx1);
            TTmat(jx,jy) = TT(kk);
        end
    end
end

surf(xmat,ymat,TTmat','facecolor','interp')
xlabel('x')
ylabel('y')
zlabel('Temperatura [^oC]')
title('Stazionario 2D')
ax=gca;
ax.FontSize=16;