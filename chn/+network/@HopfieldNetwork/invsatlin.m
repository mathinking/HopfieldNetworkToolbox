function y = invsatlin(x,u0)

    if isa(x,'gpuArray')
        y = zeros(size(x), 'gpuArray');
    else
        y = zeros(size(x));
    end
    
    y(x <= 0) = -u0;
    y(x >= 1) = u0;
    y(0 < x & x < 1) = u0 .* (2.*x(0 < x & x < 1) - 1);

end
