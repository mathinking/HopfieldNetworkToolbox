function opts = hopfieldnetOptions(varargin)
% hopfieldnetOptions Options for the HopfieldNetworkGQKP class. 
%
%   options = hopfieldnetOptions('PARAM1', VAL1, 'PARAM2', VAL2, ...) 
%   allows to define options for Generic Hopfield Networks. Then it can be 
%   used as optional options to instantiate HopfieldNetworkGQKP objects 
%   through the hopfieldnet function.
%
%   HopfieldNetworkGQKPOptions properties:
%       Setting                    - Settings used in the Simulation of
%                                    the Hopfield.
%       SimFcn                     - Simulation Function.
%       TrainFcn                   - Training Function.
%       TrainingParam              - Training Parameters.
%
%   Example:
%       Create a set of options for hopfieldnet to solve a Generic
%       Hopfield Network Problem, where the simulation function is
%       using euler's method, maximum number of iterations is 5000 and
%       the transfer function is a hyperbolic tangent transfer
%       function.
%       
%       options = hopfieldnetOptions('SimFcn', 'euler', ...
%                                    'MaxIter', 5000, ...
%                                    'TransferFcn', 'tanh');
% 
%   See also hopfieldnet, options.HopfieldNetworkGQKPOptions.

    opts = options.HopfieldNetworkGQKPOptions(varargin{:});

end
