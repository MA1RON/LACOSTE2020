function ff = L2E2(xx)

ff = [xx(1,:) + xx(2,:).^2; 4*xx(1,:) - 1./xx(2,:)];

end

