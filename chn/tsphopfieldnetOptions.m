function opts = tsphopfieldnetOptions(varargin)
% tsphopfieldnetOptions Options for the HopfieldNetworkTSP class. 
%
%   options = tsphopfieldnetOptions('PARAM1', VAL1, 'PARAM2', VAL2, ...) 
%   allows to define options for Hopfield Networks to solve the 
%   Traveling Salesman Problem (TSP). Then it can be used as optional 
%   options to instantiate HopfieldNetworkTSP objects through the 
%   tsphopfieldnet function.
%
%   HopfieldNetworkTSPOptions properties:
%       Cities                     - Defition of the TSP problem.
%       Setting                    - Settings used in the Simulation of
%                                    the Hopfield.
%       SimFcn                     - Simulation Function.
%       TrainFcn                   - Training Function.
%       TrainingParam              - Training Parameters.
%
%   Example:
%       Create a set of options for tsphopfieldnet to solves a TSP
%       with 10 cities, where city coordinates are chosen randomly,
%       distance is euclidean, maximum number of iterations is 1000, 
%       the execution environment is the cpu and the transfer function 
%       is a saturating linear transfer function.
%
%       options = tsphopfieldnetOptions('Coordinates', rand(10,2), ...
%                                       'DistanceType', 'euc', ...
%                                       'MaxIter', 1000, ...
%                                       'ExecutionEnvironment', 'cpu', ...
%                                       'TransferFcn', 'satlin');
% 
%   See also tsphopfieldnet, options.HopfieldNetworkTSPOptions.

    opts = options.HopfieldNetworkTSPOptions(varargin{:});

end
