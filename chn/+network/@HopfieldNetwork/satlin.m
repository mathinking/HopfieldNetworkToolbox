function y = satlin(x,u0)

    if isa(x,'gpuArray')
        y = zeros(size(x), 'gpuArray');
    else
        y = zeros(size(x));
    end
    y(x>-u0 & x<u0) = 0.5 * (x(x>-u0 & x<u0)./u0 + 1);
    y(x>=u0) = 1;

end
