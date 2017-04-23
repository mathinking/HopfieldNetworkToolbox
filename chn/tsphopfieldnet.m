function opts = tsphopfieldnet(networkSize, networkParameter, options)
% tsphopfieldnet Continuous Hopfield Network (CHN) to solve the TSP
%
%   net = tsphopfieldnet(networkSize, networkParameter) creates a Hopfield
%   Network to solve the TSP with networkSize cities (located in the
%   vertices of a regular polygon of size 1) using the free parameter
%   given by networkParameter.
% 
%   net = tsphopfieldnet(networkSize, networkParameter, options) specifies 
%   optional parameters that allow to configure the Hopfield Network.
%
%   tsphopfieldnet input arguments:
%       networkSize                - Number of cities in the TSP 
%                                    problem.
%       networkParameter           - Free parameter (C) in the Hopfield
%                                    Network when applied to the TSP.
%       options                    - (optional) Options for simulation 
%                                    of the Hopfield Network (object 
%                                    from the HopfieldNetworkTSPOptions
%                                    class). This options object can be
%                                    created using the
%                                    tsphopfieldnetOptions function.
% 
%   Example:
%       Build a Hopfield Network to solve a TSP with 10 cities using a 
%       value of the free parameter equal to 1e-5.
%
%       net = tsphopfieldnet(10, 1e-5);
%
%   See also tsphopfieldnetOptions, options.HopfieldNetworkTSPOptions, 
%   network.HopfieldNetworkTSP.

    if nargin == 2
        opts = network.HopfieldNetworkTSP(networkSize, networkParameter);
    elseif nargin == 3
        opts = network.HopfieldNetworkTSP(networkSize, networkParameter, options);
    else
        error('HopfieldNetworkTSP:IncorrectInputArguments','Please provide proper input arguments to HopfieldNetworkTSP: networkSize, C (and options)');
    end
    
end
