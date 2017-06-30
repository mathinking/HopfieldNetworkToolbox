classdef HopfieldNetworkTSP < network.HopfieldNetwork
    % HopfieldNetworkTSP Continuous Hopfield Network (CHN) to solve the TSP
    %
    %   This class is used to create Hopfield Networks that can solve the 
    %   Traveling Salesman Problem (TSP).
    %
    %   HopfieldNetworkTSP input arguments:
    %       networkSize                - Number of cities in the TSP 
    %                                    problem.
    %       networkParameter           - Free parameter (C) in the Hopfield
    %                                    Network when applied to the TSP.
    %       options                    - Options for simulation of the
    %                                    Hopfield Network (object from the
    %                                    tsphopfieldnetOptions class).
    % 
    %   Example:
    %       Build a Hopfield Network to solve a TSP with 10 cities using a 
    %       value of the free parameter equal to 1e-5.
    %
    %       net = network.HopfieldNetworkTSP(10, 1e-5);
    %
    %   See also tsphopfieldnetOptions, tsphopfieldnet.
    
    properties (GetAccess = private, SetAccess = private)
        Name = 'HopfieldNetworkTSP';
        Cities;
        Scheme;
    end
    
    methods (Sealed = true)
        function net = HopfieldNetworkTSP(networkSize, networkParameter, opts)

            if nargin < 2 || nargin > 3
                error('HopfieldNetworkTSP:IncorrectInputArguments','Please provide proper input arguments to HopfieldNetworkTSP: networkSize, C (and options)');
            else
                % networkSize assertion is done in hopfieldnetwork constructor. networkParameter assertion here:                                
                assert(isa(networkParameter,'float') && isreal(networkParameter) && isnumeric(networkParameter) && isscalar(networkParameter) && networkParameter > 0, ...
                    'HopfieldNetworkTSP:invalid_value', 'networkParameter (C) must be a floating point scalar greater than 0.');
            end
            
            if nargin < 3
                opts = options.HopfieldNetworkTSPOptions(); % Default Options
            else
                assert(isa(opts, 'options.HopfieldNetworkTSPOptions'), 'HopfieldNetworkTSP:invalid_datatype',...
                    'Input ''options'' must be of the class ''HopfieldNetworkTSPOptions''');
            end
            
            % Creating net object from parent class           
            net = net@network.HopfieldNetwork(networkSize,opts);
            
            net = setOptions(net, opts);
                      
            net.TrainParam.C = networkParameter;
                
            % Initialize results;
            net = init(net);
        end
        
    end
    
    % --- Methods definitions --- %    
    
	methods (Hidden = true, Access = private)
        net = setOptions(net,options);
        [net,V,U] = fixedCities(net);
        thisCity = city(net,cityPosition);
        [modifiedDistance,Ng] = neighbourDistance(net, tau_or_p);
    end

	methods (Hidden = true, Access = protected)
        net = init(net);
    end
    
    methods (Hidden = true, Access = public)
        net = reinit(net);
    end        
       
	methods (Static = true, Hidden = true, Access = private)
        p = modulo(m,n);
        loggingV(checkpointFilename,maxIter,iter,N,V,dU);
        chain = createChain(cities);
        verifyIfValidSubtours(subtours, startingPositions, myCities);
    end
    
    methods (Static = true, Hidden = true, Access = protected)

    end
    
    methods (Static = true, Hidden = false, Access = public)
        coords = polygonCoords(l,N);        
%         options = createOptions(varargin);
        myText = cityTextGeneration(N);
        d = computeDistance(coords, type);
        storeResultsinPDF(results,simulations,N_or_Problem,C,resultsPath)
    end
    
	methods(Access = public) % --- Get-Set methods --- %     
        name = getName(net);
        cities = getCities(net,field);
        setCities(net,property,value); 
        setScheme(net,scheme);        
    end
end
