clear all;close all;clc
%domain
Lx1=0.1;
Lx2=0.2;
Ly2=0.1;
Ly1=0.05;

%BC
%Tsouth=0;
Tnorth=5;
Twest=5;
Tin=20;
%conductivity
cond=10;
cp=450;
rho=7800;


%build the grid
nx2=21;nx1=11;
ny2=21;ny1=11;
xvec=linspace(0,Lx2,nx2);
yvec=linspace(0,Ly2,ny2);

[xmat,ymat]=meshgrid(xvec,yvec);
ntot=nx1*(ny1-1)+(ny2-ny1+1)*nx2;

%%impariamo a fare il grafico di una variabile sul dominio computazionale
    for jj=1:ny1-1
        temp1(1:nx1,jj)=Tnorth;
        temp1(nx1+1:nx2,jj)=NaN;
    end
    for jj=ny1:ny2
        temp1(1:nx2,jj)=Tnorth;
    end
surf(xmat,ymat,temp1');


%return
%prepare the matrix of the coefficients and the rigt-hand side
AA=sparse([],[],[],ntot,ntot,5*ntot);
rhs=zeros(ntot,1);

%build the matrix (ipotizziamo griglia uniforme)
dx=xvec(2)-xvec(1);
dy=yvec(2)-yvec(1);

%costruisco la matrice dei coefficienti in due passi
%primo loop, sul sottodominio piccolo (in basso a sx)
for icenter=2:nx1-1
   for jcenter=2:ny1
    kcenter=(jcenter-1)*nx1+icenter;
    if jcenter<ny1
        knorth=kcenter+nx1;
    else
        knorth=kcenter+nx2;
    end
    ksouth=kcenter-nx1;
    keast=kcenter+1;
    kwest=kcenter-1;
    
    AA(kcenter,kcenter)=2*cond*(1/dx^2+1/dy^2);
    AA(kcenter,knorth)=-cond/dy^2;
    AA(kcenter,ksouth)=-cond/dy^2;
    AA(kcenter,keast)=-cond/dx^2;
    AA(kcenter,kwest)=-cond/dx^2;
  
   end
end
%secondo loop, sul sottodominio grande (in alto)
for icenter=2:nx2-1
   for jcenter=1:ny2-ny1-1
    kcenter=ny1*nx1+(nx2-nx1)+(jcenter-1)*nx2 +icenter;
    knorth=kcenter+nx2;
    ksouth=kcenter-nx2;
    keast=kcenter+1;
    kwest=kcenter-1;
    
    AA(kcenter,kcenter)=2*cond*(1/dx^2+1/dy^2);
    AA(kcenter,knorth)=-cond/dy^2;
    AA(kcenter,ksouth)=-cond/dy^2;
    AA(kcenter,keast)=-cond/dx^2;
    AA(kcenter,kwest)=-cond/dx^2;
  
   end
end

%%qui inizia il blocco in cui impongo le BC
%west border
icenter=1;
for jcenter=1:ny1-1
    kcenter=(jcenter-1)*nx1+icenter;
    AA(kcenter,kcenter)=1.0;
    rhs(kcenter)=Twest;
end
%pause
for jcenter=1:(ny2-ny1+1)
    kcenter=(ny1-1)*nx1+(jcenter-1)*nx2+icenter;
    AA(kcenter,kcenter)=1.0;
    rhs(kcenter)=Twest;
end
%pause
%in border orizontal
jcenter=ny1;
for icenter=1:nx2-nx1;
    kcenter=jcenter*nx1+icenter;
    AA(kcenter,kcenter)=1.0;
    rhs(kcenter)=Tin;
end
%pause
%in border vertical
icenter=nx1;
for jcenter=1:ny1
    kcenter=(jcenter-1)*nx1+icenter;
    AA(kcenter,kcenter)=1.0;
    rhs(kcenter)=Tin;
end
%pause
%east (symmetry) border
icenter=nx2;
for jcenter=2:ny2-ny1
    kcenter=ny1*nx1+(nx2-nx1)+(jcenter-1)*icenter;
    knorth=kcenter+nx2;
    kwest=kcenter-1;
    ksouth=kcenter-nx2;
    AA(kcenter,kcenter)=2*cond/dy^2+2*cond/dx^2;
    AA(kcenter,knorth)=-cond/dy^2;
    AA(kcenter,ksouth)=-cond/dy^2;
    AA(kcenter,kwest)=-2*cond/dx^2;
    rhs(kcenter)=0;
end
%pause
% south(symmetry) border:
jcenter=1;
for icenter=2:nx1-1
    kcenter=(jcenter-1)*nx1+icenter;
    knorth=kcenter+nx1;
    keast=kcenter+1;
    kwest=kcenter-1;    
    AA(kcenter,kcenter)=2*cond/dx^2+2*cond/dy^2;
    AA(kcenter,knorth)=-2*cond/dy^2;
    AA(kcenter,keast)=-cond/dx^2;
    AA(kcenter,kwest)=-cond/dx^2;
    rhs(kcenter)=0;
end
%pause
% north:
jcenter=ny2;
for icenter=2:nx2
    kcenter=(ny1-1)*(nx1)+(jcenter-ny1)*nx2+icenter;
    AA(kcenter,kcenter)=1.0;
    rhs(kcenter)=Tnorth;
end

%%qui sotto controllo di avere imposto bene le BC
rr1=reshape(rhs(1:nx1*(ny1-1)),nx1,ny1-1);
rr2=reshape(rhs(nx1*(ny1-1)+1:end),nx2,ny2-ny1+1);
rrr=[[rr1;NaN*ones(nx2-nx1,ny2-ny1)],rr2];
surf(xmat,ymat,rrr')

%%risolvo il problmema e plotto la soluzione
TT=AA\rhs;

TTmatrix1=reshape(TT(1:nx1*(ny1-1)),nx1,ny1-1);
TTmatrix2=reshape(TT(nx1*(ny1-1)+1:end),nx2,ny2-ny1+1);
TTmatrix=[[TTmatrix1;NaN*ones(nx2-nx1,ny2-ny1)],TTmatrix2];

figure
surf(xmat,ymat,TTmatrix','facecolor','interp')
xlabel('x','fontsize',16)
ylabel('y','fontsize',16)
zlabel('Temperature','fontsize',16)
%%fine del conto stazionario

%%calcolo il flusso (integro sui due bordi W e N)
fluxleft=0;
for i=1:ny1
    if i==1;
        delta=dy/2;
    else
        delta=dy;
    end
    fluxleft=fluxleft+cond/dx*(TT((i-1)*nx1+2)-TT((i-1)*nx1+1))*delta;
end
for i=1:(ny2-ny1)
    if i==(ny2-ny1);
        delta=dy/2;
    else
        delta=dy;
    end
    fluxleft=fluxleft+cond/dx*(TT((ny1-1)*nx1+i*nx2+2)-TT((ny1-1)*nx1+i*nx2+1))*delta;
end
fluxtop=0;
nn=0;
for i=nx1*(ny1-1)+nx2*(ny2-ny1)+1:length(TT)
    nn=nn+1;
    if nn==1 | nn==nx2
        delta=dx/2;
    else
        delta=dx;
    end
    fluxtop=fluxtop+cond/dy*(TT(i)-TT(i-(2*nn+1)))*delta;
end
fluxtot=fluxleft+fluxtop

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%transient part
timescale=rho*cp/cond;
deltat=1800;
BB=diag(ones(length(AA),1),0)+AA/rho/cp*deltat;
%BB=AA;
rhs=zeros(length(BB),1);

%%blocco in cui assegno le condizoni al contorno su tutti i bordi nella
%modificando SOLO gli elementi della matrice BB (per l'equivalente nel
%vettore dei termini noti, devo aspettare di essere nel loop che avanza in
%tempo
%west border
icenter=1;
for jcenter=1:ny1-1
    kcenter=(jcenter-1)*nx1+icenter;
    BB(kcenter,kcenter)=1.0;
end
%pause
for jcenter=1:(ny2-ny1+1)
    kcenter=(ny1-1)*nx1+(jcenter-1)*nx2+icenter;
    BB(kcenter,kcenter)=1.0;
end
%pause
%in border orizontal
jcenter=ny1;
for icenter=1:nx2-nx1;
    kcenter=jcenter*nx1+icenter;
    BB(kcenter,kcenter)=1.0;
end
%pause
%in border vertical
icenter=nx1;
for jcenter=1:ny1
    kcenter=(jcenter-1)*nx1+icenter;
    knorth=kcenter+nx1;
    kwest=kcenter-1;
    ksouth=kcenter-nx1;
    BB(kcenter,kcenter)=1;
end
%pause
%east (symmetry) border
icenter=nx2;
for jcenter=2:ny2-ny1
    kcenter=ny1*nx1+(nx2-nx1)+(jcenter-1)*icenter;
    knorth=kcenter+nx2;
    kwest=kcenter-1;
    ksouth=kcenter-nx2;
    BB(kcenter,kcenter)=2*cond/dy^2+2*cond/dx^2;
    BB(kcenter,knorth)=-cond/dy^2;
    BB(kcenter,ksouth)=-cond/dy^2;
    BB(kcenter,kwest)=-2*cond/dx^2;
end
% south(symmetry) border:
jcenter=1;
for icenter=2:nx1-1
    kcenter=(jcenter-1)*nx1+icenter;
    knorth=kcenter+nx1;
    keast=kcenter+1;
    kwest=kcenter-1;    
    BB(kcenter,kcenter)=2*cond/dx^2+2*cond/dy^2;
    BB(kcenter,knorth)=-2*cond/dy^2;
    BB(kcenter,keast)=-cond/dx^2;
    BB(kcenter,kwest)=-cond/dx^2;
end
%pause
% north:
jcenter=ny2;
for icenter=2:nx2
    kcenter=(ny1-1)*(nx1)+(jcenter-ny1)*nx2+icenter;
    BB(kcenter,kcenter)=1.0;
end
%fine del blocco per assgenare le BC

%%da qui preparo le variabili per il loop in tempo
TTold=TT;
qv=zeros(size(TTold));
time=0;
Tmid=TT(nx1*(ny1-1)+ceil(nx1/2));
driver=Tnorth;
figure
xlabel('x','fontsize',16)
ylabel('y','fontsize',16)
zlabel('Temperature (^oC)','fontsize',16)
zlim([-10 20])

%%per il video (sfizio non richiesto)
h.visible='off';
loops = 42;
v = VideoWriter('ex1bp.avi');
open(v);
M(loops) = struct('cdata',[],'colormap',[]);

%%loop in tempo (uso ciclo for, perchè risolvo per un periodo di 24h 
%quindi so a priori quante volte entrare nel ciclo)
for ii=2:round(24*3600/deltat)  %il tempo nella formula era in ore
    time(ii)=time(ii-1)+deltat;  %calcolo l'istante di tempo in cui mi trovo
    Tnorth=5+15*sin(2*pi*time(ii)/3600/24); % calcolo la BC a quel tempo...
    driver(ii)=Tnorth;  %...e la assegno alla variabile (vettore) driver
    Twest=Tnorth;
    rhs=(TTold+qv*deltat); %incomincio a creare il vettore dei termini noti
%%da qui il blocco per inserire al termine noto le BC
    %west border
        icenter=1;
        for jcenter=1:ny1-1
        kcenter=(jcenter-1)*nx1+icenter;
        rhs(kcenter)=Twest;
        end
        for jcenter=1:(ny2-ny1+1)
        kcenter=(ny1-1)*nx1+(jcenter-1)*nx2+icenter;
        rhs(kcenter)=Twest;
        end
    %in border orizontal
        jcenter=ny1;
        for icenter=1:nx2-nx1;
        kcenter=jcenter*nx1+icenter;
        rhs(kcenter)=Tin;
        end
    %pause
    %in border vertical
        icenter=nx1;
        for jcenter=1:ny1
        kcenter=(jcenter-1)*nx1+icenter;
        rhs(kcenter)=Tin;
        end
    %pause
    %east (symmetry) border
        icenter=nx2;
        for jcenter=2:ny2-ny1
        kcenter=ny1*nx1+(nx2-nx1)+(jcenter-1)*icenter;
        rhs(kcenter)=0;
        end
    % south(symmetry) border:
        jcenter=1;
        for icenter=2:nx1-1
        kcenter=(jcenter-1)*nx1+icenter;
        rhs(kcenter)=0;
    end
    %pause
    % north:
        jcenter=ny2;
        for icenter=2:nx2
        kcenter=(ny1-1)*(nx1)+(jcenter-ny1)*nx2+icenter;
        rhs(kcenter)=Tnorth;
    end
%%fine del blocco per inserire le BC

%%risolvo il rpoblema
    TTnew=BB\rhs;
    
%%avanzo in tempo    
    TTold=TTnew;

%calcolo la temperatura nel punto di interesse 
%mi sono stufata e determino gli indici alla spiccia
    Tmid(ii)=TTnew(nx1*(ny1-1)+ceil(nx1/2));
    
%faccio il reshaper pe poter plottar ela soluzione
%inserisco i NaN dove servono
    TTmatrix1=reshape(TTnew(1:nx1*(ny1-1)),nx1,ny1-1);
    TTmatrix2=reshape(TTnew(nx1*(ny1-1)+1:end),nx2,ny2-ny1+1);
    TTmatrix=[[TTmatrix1;NaN*ones(nx2-nx1,ny2-ny1)],TTmatrix2];

    
    %figure
    surf(xmat,ymat,TTmatrix','facecolor','interp')
    xlabel('x (m)')
    ylabel('y (m)')
    zlabel('Temperature (^oC)')

    hold on
    zlim([-10 20])
%questo per fare il video dell'evoluzione della mappa di T nel tempo (non
%richiesto)
    ax = gca;
    ax.NextPlot = 'replaceChildren';
    drawnow
    M(ii-1) = getframe;
%     pause
end
%questo per fare il video dell'evoluzione della mappa di T nel tempo (non
%richiesto)
h.Visible = 'on';
% movie(M)
writeVideo(v,M)
close(v)

%la figura sotto per la diagnostica richiesta
 figure
 plot(time/3600,driver,'k--','Linewidth',2)
 hold on
 plot(time/3600,Tmid,'r-','Linewidth',2)
 xlabel('Time (h)')
 ylabel('Temperature (^oC)')
 legend('T_o_u_t','T(0.5 m 0.5 m)')
 grid on