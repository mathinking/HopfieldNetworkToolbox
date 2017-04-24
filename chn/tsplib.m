function problem = tsplib(name)
% tsplib Define a TSPLIB problem
%
%   problem = tsplib(name) creates a TSPLIB [1] problem with relevant
%   information to solve a TSP instance.
% 
%   [1] Reinelt, G. (1991). TSPLIB—A traveling salesman problem library. 
%       ORSA journal on computing, 3(4), 376-384.
%
%   tsplib input arguments:
%
%       name                       - TSPLIB problem name.
%
%   Example:
%       Create a TSPLIB class for the 'berlin52' TSPLIB problem:
% 
%       problem = tsplib({'berlin52'});
%
%   See also tsphopfieldnetOptions, tsphopfieldnet.

    problem = utils.TSPLIB(name);

end
