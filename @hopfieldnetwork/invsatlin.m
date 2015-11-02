function y = invsatlin(x,u0)
    y = zeros(size(x));
    for j = 1:size(x,2)
        for i = 1:size(x,1);
            if x(i,j) < 0
                y(i,j) = u0;
            elseif x(i,j) > 1
                y(i,j) = -u0;
            else
                y(i,j) = u0 * (2*x(i,j) - 1);
            end       
        end
    end
end
