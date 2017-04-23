function d = computeDistance(coords, type)
% COMPUTEDISTANCE obtains the matrix of distances for all points of the
% COORDS matrix
%
%   d = COMPUTEDISTANCE(COODS, TYPE) computes the square and symmetric matrix of
%   distances D (using the Euclidean distance) for all coordinate points in
%   the matrix COORDS.
% 
%   See also PDIST, SQUAREFORM

    switch lower(type)
        case 'geo'
            d = cityDistances_GEO(coords);
        case 'euc_2d'
            d = cityDistances_EUC_2D(coords); % Uses round (for TSPLIB)
        case 'euc'
            d = pdist(coords,'euclidean');
            d = squareform(d);        
        case 'att'
            d = cityDistances_ATT(coords);
        case 'ceil_2d'
            d = cityDistances_CEIL_2D(coords);
        otherwise
            error('Unknown distance matrix case');
    end
    
end

function d = cityDistances_GEO(coords)

    N = size(coords,1);
    d = zeros(N);
    RRR = 6378.388;
    latitude = coords(:,1); 
    longitude = coords(:,2);
    for i = 1:N
        for j = (i+1):N
            q1 = cos( longitude(i) - longitude(j) );
            q2 = cos( latitude(i) - latitude(j) );
            q3 = cos( latitude(i) + latitude(j) ); 
            d(i,j) =...
                fix( RRR * acos( 0.5*((1+q1)*q2 - (1-q1)*q3) ) + 1);
            d(j,i) = d(i,j);
        end
    end
    
end

function d = cityDistances_EUC_2D(coords)

    d = pdist(coords,'euclidean');
    d = round(squareform(d));
    
end

function d = cityDistances_ATT(coords)

    N = length(coords);
    d = zeros(N);

    for i = 1:N
        for j = (i+1):N
            xd = coords(i,1) - coords(j,1);
            yd = coords(i,2) - coords(j,2);
            r = sqrt( (xd.^2 + yd.^2)./10 );
            t = round(r);

            if t < r
                d(i,j) = t + 1;
            else
                d(i,j) = t;
            end
            d(j,i) = d(i,j);        
        end
    end    

end

function d = cityDistances_CEIL_2D(coords)

    d = pdist(coords,'euclidean');
    d = ceil(squareform(d));

end
