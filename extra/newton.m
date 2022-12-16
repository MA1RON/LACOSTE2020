function [xm, fm, jj] = newton(x0, ff, tol, jmax, showplot)
    accur = 1.1;
    dx = ff(x0)/accur;
    xx = linspace(x0-dx*accur^2,x0+dx*accur^2);
    df = @(x) (ff(x+dx)-ff(x-dx))/2/dx;
    
    if showplot
        plot(xx, ff(xx))
        hold on
        plot(xx, df(xx))
        hold on
        plot([xx(1), xx(end)], [0, 0], 'k-.')
    end
    err = 1;
    jj = 1;
    xm = x0;
    while err>tol
        xm = xm - ff(xm)/df(xm);
        
        if showplot
            pause(.3)
            plot(xm,ff(xm),'ko', 'markerfacecolor','y','markersize',8)
            hold on
        end
        
        if jj == jmax
            warning('la soluzione non converge')
            return
        end
        
        err = norm(ff(xm));
    end
end