function jac = jacbuilder(TT,kk,kder,kder2,hh,qqq,dx)
nn = length(qqq);

ki = kk(TT(1));
kp = kk(TT(2:end-1));
ks = kk(TT(end));

dki = kder(TT(1));
dkp = kder(TT(2:end-1));
dks = kder(TT(end));

dkp2 = kder2(TT(2:end-1));

diag_s = [0;...
          ki/hh/dx;... % coco
          dkp.*(TT(3:end)-TT(1:end-2))/2/dx^2 + kp/dx^2];
diag_p = [dki/hh/dx*(TT(2)-TT(1))-ki/hh/dx-1;... % coco
          dkp2.*(TT(3:end)-TT(1:end-2)).^2/4/dx^2 + dkp.*(TT(3:end)-2*TT(2:end-1)+TT(1:end-2))/dx^2 + kp.*(-2/dx^2);...
          dks/hh/dx*(TT(end-1)-TT(end))-ks/hh/dx-1]; % coco
% diag_p = [-ki/hh/dx-1;... % coco
%           dkp2.*(TT(3:end)-TT(1:end-2)).^2/4/dx^2 + dkp.*(TT(3:end)-2*TT(2:end-1)+TT(1:end-2))/dx^2 + kp.*(-2/dx^2);...
%           -ks/hh/dx-1]; % coco
diag_i = [-dkp.*(TT(3:end)-TT(1:end-2))/2/dx^2 + kp/dx^2;...
          ks/hh/dx;... % coco
          0];

jac = spdiags([diag_i diag_p diag_s],-1:1,nn,nn);
end