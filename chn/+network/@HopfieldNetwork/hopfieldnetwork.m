classdef HopfieldNetwork < handle
    %HopfieldNetwork Continuous Hopfield Network (CHN)
    %
    %   This class defines abstract methods and properties for Hopfield 
    %   Networks. Therefore it can only be instantiated from within a child 
    %   class.
    %
    %   Input arguments properties:
    %       Setting                    - Settings used in the Simulation of
    %                                    the Hopfield Network.
    %       SimFcn                     - Simulation Function.
    %       TrainFcn                   - Training Function.
    %       TrainingParam              - Training Parameters.
    %       Results                    - Simulation Results.
    %
    %   See also options.HopfieldNetworkOptions, hopfieldnet, tsphopfieldnet.
    
    properties (GetAccess = protected, SetAccess = protected)
        Setting;
        SimFcn;
        TrainFcn;
        TrainParam;
        Results;
    end
    
    methods (Sealed = true)
        function net = HopfieldNetwork(networkSize, opts)
            
            % Input checking
            if nargin < 1 || nargin > 2
                error('HopfieldNetwork:IncorrectInputArguments','Please provide proper input arguments to HopfieldNetwork: networkSize (and options)');
            else
                assert(isreal(networkSize) && isnumeric(networkSize) && isscalar(networkSize) && all(mod(networkSize,1)==0) && networkSize >= 1,...
                    'HopfieldNetwork:invalid_value', [num2str(networkSize), ' must be a scalar integer greater or equal than 1.']);
                if nargin == 2
                    assert(isa(opts,'options.HopfieldNetworkOptions') || isa(opts,'options.HopfieldNetworkGQKPOptions') || isa(opts,'options.HopfieldNetworkTSPOptions'),...
                        'HopfieldNetwork:invalid_datatype', 'options must be an object from the class HopfieldNetworkOptions, HopfieldNetworkGQKPOptions or HopfieldNetworkTSPOptions.');
                end
            end
            
            if ~isa(networkSize,'float')
                networkSize = double(networkSize);
            end
                
            % Setting network attributes
            net.TrainParam.N = networkSize;

            % Setting default options if not provided
            if nargin < 2 
                opts = options.HopfieldNetworkOptions();               
            end      

            % Setting desired options available in options structure
            net = setOptions(net, opts);
            
            net.Setting = orderfields(net.Setting);

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
        setTrainParam(net,property,value);        
    end
        
	methods (Hidden = true, Access = private)
        %net = addDefaultOptionValues(net, options);
        net = setOptions(net, options);
    end

    methods (Hidden = false, Access = protected)
        net = init(net);
    end
    
%     methods (Static = true, Access = public)
%         options = createOptions(varargin);
%     end
    
    methods (Static = true, Access = protected)
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
        setResults(net,property,value);        
    end
end
