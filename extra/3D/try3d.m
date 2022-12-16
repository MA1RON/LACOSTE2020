nn = 11;

Lx = 1;
Ly = 2;
Lz = 1;

xvet = linspace(0,Lx,nn);
yvet = linspace(0,Ly,nn);
zvet = linspace(0,Lz,nn);

t = ones(1,nn^3);

for jx = 1:nn
    for jy = 1:nn
        for jz = 1:nn
            kk = jx + (jy-1)*nn + (jz-1)*nn^2;
            scatter3(jx,jy,jz,20,t(kk));
            hold on
        end
    end
end
colorbar