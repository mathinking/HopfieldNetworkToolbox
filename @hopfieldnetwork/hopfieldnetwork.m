classdef hopfieldnetwork < handle
    %HOFIELDNETWORK Continuous Hopfield Network (CHN)
    
    properties (GetAccess = protected, SetAccess = protected)
        trainFcn;
        trainParam;
        setting;
        simFcn;
%         id; % If required, do get-set functions
        results;
    end
    
    methods (Sealed = true)
        function net = hopfieldnetwork(networkSize, options)
            
            % Input checking
            if nargin < 1 || nargin > 2
                error('hopfieldnetwork:IncorrectInputArguments','Please provide proper input arguments to hopfieldnetwork: networkSize (and options)');
            else
                assert(isa(networkSize,'double'), 'hopfieldnetwork:invalid_datatype', [networkSize, ' must be double.']);
                assert(all(networkSize >= 1), 'hopfieldnetwork:invalid_value', '''networkSize'' must be greater or equal to 1.');
                if nargin == 2
                    assert(isstruct(options), 'hopfieldnetwork:invalid_datatype', 'options must be a structure created using createOptions.');
                end
            end

            % Setting network attributes
            net.trainParam.N = networkSize;

            % Setting default options if not provided
            if nargin < 2 
                options = hopfieldnetwork.createOptions();
                net = addDefaultOptionValues(net, options);

            % Setting desired options available in options structure
            else    
                net = setOptions(net, options);
                net = addDefaultOptionValues(net, options);
            end
            
            net.setting = orderfields(net.setting);

        end
    end
    
    % --- Methods definitions --- %    
    methods (Abstract = true)
        net = train(net);
        net = sim(net);
        disp(net);
        plot(net, varargin);
        
        % --- Set methods --- %
        setTrainFcn(net,trainFcn);
        setSimFcn(net,simFcn);
    end
        
	methods (Hidden = true, Access = private)
        net = addDefaultOptionValues(net, options);
        net = setOptions(net, options);
    end

    methods (Hidden = false, Access = protected)
        net = init(net);
    end
    
    methods (Static = true, Access = public)
        options = createOptions(varargin);
    end
    
    methods (Static = true, Access = protected)
        isValidSetting(property,value);
        y = satlin(x,u0);
        y = invsatlin(x,u0);        
    end 
    
    methods(Access = public) % --- Get-Set methods --- %
        trainFcn = getTrainFcn(net);
        trainParam = getTrainParam(net,field);
        setting = getSetting(net,field);
        simFcn = getSimFcn(net);
        results = getResults(net,field);
        setSetting(net,property,value);
        setTrainParam(net,property,value);
        setResults(net,property,value);        
    end
end
