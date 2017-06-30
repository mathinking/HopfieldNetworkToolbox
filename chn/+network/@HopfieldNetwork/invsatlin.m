function y = invsatlin(x,u0)
    y = zeros(size(x));
    
    y(x <= 0) = -u0;
    y(x >= 1) = u0;
    y(0 < x & x < 1) = u0 .* (2.*x(0 < x & x < 1) - 1);

end
