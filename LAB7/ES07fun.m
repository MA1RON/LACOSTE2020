function [ff] = ES07fun(tt,yy,alpha,beta,gamma,delta)

    ff = zeros(2,1);
    ff(1) = alpha*yy(1)-beta*yy(1)*yy(2);
    ff(2) = gamma*yy(1)*yy(2)-delta*yy(2);


end

