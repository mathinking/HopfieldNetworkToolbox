function coords = polygonCoords(l,N)
% POLYGONCOORDS computes the coordinates of the vertices of a polygon with 
% N edges. 
%
%   COORDS = POLYGONCOORDS(L,N) calculates a two column vector COORDS with
%   the X and Y coordinates of a polygon with N edges, each of length L 
%   from the center of the POLYGON to its edges.

    t = (0:pi/(N/2):2*pi-pi/(N/2))';
    coords = [l*sin(t),l*cos(t)];    
end
