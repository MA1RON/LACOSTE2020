%% LaCoSTe 2020/2021 
%% Lab 9 Ex 5 - Soluzioni
%% (variazione dimensione pozzo)

clear all
close all
clc

%Dimensioni
ws=2e-1;
wf=2e-1;
%Proprietà e temperature
cond=400;
T_h=25;
h=3e4;
T_c=75;

%% Costruzione dominio computazionale
% Nomino i vertici del dominio A,B,C,D,E,F partendo dal vertice in basso a
% sinistra e procedendo in senso antiorario.
% F---E
% |   |
% | C-D
% | |
% A-B

avec=3e-2:1e-2:20e-2; % Per variazione dimensioni pozzo (terza richiesta)

%% Iterazioni con diverse dimensioni del pozzo
for i=1:length(avec)
    dx=1e-2;
    dy=dx;
    
    a=avec(i); %Calcolo dimensioni a e b per le diverse iterazioni
    b=4e-1-a;

    xvec=(0:dx:(ws/2+wf/2))';
    yvec=(0:dy:(a+b))';
    [xmat,ymat]=meshgrid(xvec,yvec);

    %Calcolo numero di nodi sui diversi bordi
    n1=length(0:dx:ws/2); %AB
    n2=length(xvec); %EF
    n3=length(yvec); %FA
    n4=length(0:dy:b); %BC

    nvec(i)=n1; %Mi serve per lo studio di grid indipendence

    ntot=n2*n3-(n2-n1)*(n4-1);
    A=sparse([],[],[],ntot,ntot);
    rhs=zeros(ntot,1);

    %% Matrice dei coefficienti
    %Per i punti interni del rettangolo piccolo (in basso)
    for icenter=2:n1-1
       for jcenter=2:n4
        kcenter=(jcenter-1)*n1+icenter;
        if jcenter<n4
            knorth=kcenter+n1;
        else
            knorth=kcenter+n2;
        end
        ksouth=kcenter-n1;
        keast=kcenter+1;
        kwest=kcenter-1;

        A(kcenter,kcenter)=2*cond*(1/dx^2+1/dy^2);
        A(kcenter,knorth)=-cond/dy^2;
        A(kcenter,ksouth)=-cond/dy^2;
        A(kcenter,keast)=-cond/dx^2;
        A(kcenter,kwest)=-cond/dx^2;
       end
    end
    %Per i punti interni del rettangolo grande (in alto)
    for icenter=2:n2-1
       for jcenter=n4+1:n3-1
        kcenter=n1*(n4-1)+(jcenter-n4)*n2+icenter;
        knorth=kcenter+n2;
        ksouth=kcenter-n2;
        keast=kcenter+1;
        kwest=kcenter-1;

        A(kcenter,kcenter)=2*cond*(1/dx^2+1/dy^2);
        A(kcenter,knorth)=-cond/dy^2;
        A(kcenter,ksouth)=-cond/dy^2;
        A(kcenter,keast)=-cond/dx^2;
        A(kcenter,kwest)=-cond/dx^2;
       end
    end

    %% Impongo le condizioni al contorno sui bordi del dominio
    %North (EF, Dirichlet) (estremi del segmento inclusi)
    jcenter=n3;
    for icenter=1:n2
        kcenter=n1*(n4-1)+(jcenter-n4)*n2+icenter;
        A(kcenter,kcenter)=1.0;
        rhs(kcenter)=T_c;
    end
    %West (FA, Neumann adiabatico) (estremi del segmento esclusi)
    icenter=1;
    for jcenter=2:n4-1
        kcenter=(jcenter-1)*n1+icenter;
        knorth=kcenter+n1;
        ksouth=kcenter-n1;
        keast=kcenter+1;
        A(kcenter,kcenter)=2*cond/dy^2+2*cond/dx^2;
        A(kcenter,knorth)=-cond/dy^2;
        A(kcenter,ksouth)=-cond/dy^2;
        A(kcenter,keast)=-2*cond/dx^2;
        rhs(kcenter)=0;
    end
    for jcenter=n4:n3-1
        kcenter=n1*(n4-1)+(jcenter-n4)*n2+icenter;
        knorth=kcenter+n2;
        ksouth=kcenter-n2;
        keast=kcenter+1;
        A(kcenter,kcenter)=2*cond/dy^2+2*cond/dx^2;
        A(kcenter,knorth)=-cond/dy^2;
        A(kcenter,ksouth)=-cond/dy^2;
        A(kcenter,keast)=-2*cond/dx^2;
        rhs(kcenter)=0;
    end
    %South (AB, Neumann adiabatico) (estremi del segmento esclusi)
    for icenter=2:n1-1
        kcenter=icenter;
        knorth=kcenter+n1;
        keast=kcenter+1;
        kwest=kcenter-1;
        A(kcenter,kcenter)=2*cond/dy^2+2*cond/dx^2;
        A(kcenter,keast)=-cond/dx^2;
        A(kcenter,kwest)=-cond/dx^2;
        A(kcenter,knorth)=-2*cond/dy^2;
        rhs(kcenter)=0;
    end
    %South (CD, Robin) (estremi del segmento esclusi)
    for icenter=n1+1:n2-1
        kcenter=n1*(n4-1)+icenter;
        knorth=kcenter+n2;
        keast=kcenter+1;
        kwest=kcenter-1;
        A(kcenter,kcenter)=2*(cond/dy^2+cond/dx^2+h/dy);
        A(kcenter,keast)=-cond/dx^2;
        A(kcenter,kwest)=-cond/dx^2;
        A(kcenter,knorth)=-2*cond/dy^2;
        rhs(kcenter)=2*h*T_h/dy;
    end
    %East (DE, Neumann adiabatico) (estremi del segmento esclusi)
    icenter=n2;
    for jcenter=n4+1:n3-1
        kcenter=n1*(n4-1)+(jcenter-n4)*n2+icenter;
        knorth=kcenter+n2;
        ksouth=kcenter-n2;
        kwest=kcenter-1;
        A(kcenter,kcenter)=2*cond/dy^2+2*cond/dx^2;
        A(kcenter,knorth)=-cond/dy^2;
        A(kcenter,ksouth)=-cond/dy^2;
        A(kcenter,kwest)=-2*cond/dx^2;
        rhs(kcenter)=0;
    end
    %East (BC, Robin) (estremi del segmento esclusi)
    icenter=n1;
    for jcenter=2:n4-1
        kcenter=(jcenter-1)*n1+icenter;
        knorth=kcenter+n1;
        ksouth=kcenter-n1;
        kwest=kcenter-1;
        A(kcenter,kcenter)=2*(cond/dy^2+cond/dx^2+h/dx);
        A(kcenter,knorth)=-cond/dy^2;
        A(kcenter,ksouth)=-cond/dy^2;
        A(kcenter,kwest)=-2*cond/dx^2;
        rhs(kcenter)=2*h*T_h/dx;
    end

    %% Impongo le condizioni al contorno ai 4 vertici mancanti

    % Punto A (Neumann adiabatico)
    kcenter=1;
    knorth=kcenter+n1;
    keast=kcenter+1;
    A(kcenter,kcenter)=cond/dy^2+cond/dx^2;
    A(kcenter,knorth)=-cond/dy^2;
    A(kcenter,keast)=-cond/dx^2;
    rhs(kcenter)=0;

    % Punto B (Robin)
    kcenter=n1;
    knorth=kcenter+n1;
    kwest=kcenter-1;
    A(kcenter,kcenter)=-(cond/dy^2+cond/dx^2+h/dx);
    A(kcenter,knorth)=cond/dy^2;
    A(kcenter,kwest)=cond/dx^2;
    rhs(kcenter)=-h/dx*T_h;

    % Punto D (Robin)
    kcenter=n1*n4+n2-n1;
    knorth=kcenter+n2;
    kwest=kcenter-1;
    A(kcenter,kcenter)=-(cond/dy^2+cond/dx^2+h/dy);
    A(kcenter,knorth)=cond/dy^2;
    A(kcenter,kwest)=cond/dx^2;
    rhs(kcenter)=-h/dy*T_h;

    % Punto C (Robin)
    kcenter=n1*n4;
    knorth=kcenter+n2;
    ksouth=kcenter-n1;
    kwest=kcenter-1;
    keast=kcenter+1;
    A(kcenter,kcenter)=3*cond/dy^2+3*cond/dx^2+h/dx+h/dy;
    A(kcenter,knorth)=-2*cond/dy^2;
    A(kcenter,ksouth)=-cond/dy^2;
    A(kcenter,kwest)=-cond/dx^2;
    A(kcenter,keast)=-2*cond/dx^2;
    rhs(kcenter)=h*T_h/dx+h*T_h/dy;

    %% Risolvo il problema
    T=A\rhs;

    %% Calcolo la dissipazione di calore
    % Il valore q_1 [W/m] corrisponde alla dissipazione di calore del singolo
    % microcanale, per unità di lunghezza nella direzione normale alla figura.
    % Il valore q_2 [W/m2] corrisponde alla dissipazione di calore considerando
    % una superficie del chip di 10*10mm.
    q_1=2*(trapz(yvec(1:n4),h*(T(n1:n1:n1*n4)-T_h))+...
           trapz(xvec((n1+1):n2),h*(T((n1*n4+1):(n1*n4+n2-n1))-T_h)));
    L=1e-2;
    q_2(i)=q_1*L*(L/(2*a)); % Moltiplico il calore per unità di lunghezza per la lunghezza del chip in direzione normale alla figura e poi "proporziono" rispetto alla direzione orizzontale.
end

%% Plot
figure(3)
plot(avec,q_2,'Linewidth',2)
grid minor
xlabel('a [m]')
ylabel('Massimo calore dissipabile [W]')