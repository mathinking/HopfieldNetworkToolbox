function y = satlin(x,u0)
    y = zeros(size(x));
    y(x>-u0 & x<u0) = 0.5 * (x(x>-u0 & x<u0)./u0 + 1);
    y(x>=u0) = 1;
end
