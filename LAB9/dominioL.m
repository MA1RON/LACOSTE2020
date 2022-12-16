clc; clear; close all;
%% dati
% topologia
Lx1 = .2;
Lx2 = .1;
Ly1 = .1;
Ly2 = .05;

% fisici
T0 = 5;
Tint = 20;

%% griglie
nx1 = 4;
nx2 = 3;
ny1 = 4;
ny2 = 3;
ntot = nx1*ny1-(nx1-nx2)*(ny1-ny2);
xvec = linspace(0,Lx1,nx1);
yvec = linspace(0,Ly1,ny1);
[xmat, ymat] = meshgrid(xvec,yvec);
TT = ones(ntot,1)*T0;

TTmat = zeros(nx1,ny1);
for jx = 1:nx1
    for jy = 1:ny1
        if jx>nx2 && jy>ny2
            TTmat(jx,jy) = nan;
        elseif jy<=nx2
            kk = jx+(jy-1)*nx1;
            TTmat(jx,jy) = TT(kk);
        else
            kk = jx+(jy-1-ny2)*nx2+nx1*ny2;
            TTmat(jx,jy) = TT(kk);
        end
    end
end

%% surfo 
surf(xmat, ymat, TTmat, 'facecolor', 'interp')