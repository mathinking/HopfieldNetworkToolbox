classdef tsphopfieldnet < hopfieldnetwork
    %TSPHOFIELDNETWORK Continuous Hopfield Network (CHN) for TSP
    
    properties (GetAccess = private, SetAccess = private)
        name = 'Continuous Hopfield Network (CHN) for TSP';
        cities;
    end
    
    methods (Sealed = true)
        function net = tsphopfieldnet(networkSize, networkParameter, options)

            if nargin < 2 || nargin > 3
                error('tsphopfieldnet:IncorrectInputArguments','Please provide proper input arguments to tsphopfieldnet: networkSize, C (and options)');
            else
                % networkSize assertion done in hopfieldnetwork constructor
                assert(networkSize > 1, 'tsphopfieldnet:invalid_value', '''networkSize'' must be greater than 2.');
                
                assert(numel(networkParameter) == 1, 'networkParameters must be 1: ''C''.')
                
                C = networkParameter(1);
                assert(isa(C,'double'), 'tsphopfieldnet:invalid_datatype', '''C'' must be double.');
                assert(C > 0, 'tsphopfieldnet:invalid_value', '''C'' must be greater than 0.');
                
            end
                        
            if nargin < 3
                options = tsphopfieldnet.createOptions();
            end
            
            % Creating net object from parent class           
            net = net@hopfieldnetwork(networkSize,options);
            
            % Setting default options if not provided
            if nargin < 3            
                net = addDefaultOptionValues(net, options);           
                
            % Setting desired options available in options structure
            else    
                net = setOptions(net, options);
                net = addDefaultOptionValues(net, options);
            end
            
            if isempty(net.cities.d)
                net.cities.d = tsphopfieldnet.computeDistance(net.cities.coords,net.cities.type);
            end

            net.cities = orderfields(net.cities);
            
            net.trainParam.C = C;
            
            if (~isempty(net.cities.coords) && size(net.cities.coords,1) ~= net.trainParam.N) || ...
                    (isempty(net.cities.coords) && size(net.cities.d,1) ~= net.trainParam.N)
                error('Mismatch between number of coordinates, distance matrix and neurons');
            end
                
            % Initialize results;
            net = init(net);
        end
        
    end
    
    % --- Methods definitions --- %    
    
	methods (Hidden = true, Access = private)
        net = addDefaultOptionValues(net, options);
        net = setOptions(net,options);
        net = verifyIfValidSubtours(net);
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
        loggingV(iter,V,dU);
        chain = createChain(cities);
    end
    
    methods (Static = true, Hidden = true, Access = protected)
        isValidSetting(property,value);        
    end
    
    methods (Static = true, Hidden = false, Access = public)
        coords = polygonCoords(l,N);        
        options = createOptions(varargin);
        myText = cityTextGeneration(N);
        d = computeDistance(coords, type);
        storeResultsinPDF(results,simulations,N_or_Problem,C,resultsPath)
    end
    
	methods(Access = public) % --- Get-Set methods --- %     
        name = getName(net);
        cities = getCities(net,field);
        setCities(net,property,value); 
    end
end