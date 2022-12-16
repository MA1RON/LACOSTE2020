function [ff] = ES06funfun(tt,yy)

    ff = zeros(2,1);
    ff(1) = -4/3*(yy(2)-3);
    ff(2) = 3/4*(yy(1)-2);

end

