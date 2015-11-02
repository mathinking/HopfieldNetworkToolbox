classdef generichopfieldnet < hopfieldnetwork
    %GENERICHOFIELDNETWORK Continuous Hopfield Network (CHN)
    
    properties (GetAccess = public, SetAccess = public)
        name = 'Continuous Hopfield Network (CHN)';
        originalParameters;
        problemParameters;
    end
    
    methods
        function net = generichopfieldnet(P, q, A, b, Aeq, beq, options)
            % Solving the GQKP
            % min{1/2 * v' * P * v + q' * v}
            % subject to
            %
            % R * v = b
            % v_i belongs to {0,1} for all i = 1,...n
            % v_{n+k} belongs to [0,1], k = 1,..., m1
            
            if nargin < 6 || nargin > 7
                error('generichopfieldnet:IncorrectInputArguments','Please provide proper input arguments to generichopfieldnet: P, q, R, b (and options)');
            else                
                % asserting inputs
                assert(isa(P,'double'), 'generichopfieldnet:invalid_datatype', 'P must be double.');
                assert(isa(q,'double'), 'generichopfieldnet:invalid_datatype', 'q must be double.');
                assert(isa(A,'double'), 'generichopfieldnet:invalid_datatype', 'A must be double.');
                assert(isa(b,'double'), 'generichopfieldnet:invalid_datatype', 'b must be double.');
                assert(isa(Aeq,'double'), 'generichopfieldnet:invalid_datatype', 'Aeq must be double.');
                assert(isa(beq,'double'), 'generichopfieldnet:invalid_datatype', 'beq must be double.');                
            end           
                        
            if nargin < 7
                options = generichopfieldnet.createOptions();
            end
            
            % Creating net object from parent class
            networkSize = size(P,1) + size(A,1);
            net = net@hopfieldnetwork(networkSize,options);

            net.originalParameters.P = P;
            net.originalParameters.q = q;
            net.originalParameters.A = A;
            net.originalParameters.b = b;
            net.originalParameters.Aeq = Aeq;
            net.originalParameters.beq = beq;
                       
            net.problemParameters.n  = size(net.originalParameters.P,1);
            net.problemParameters.m1 = size(net.originalParameters.A,1);
            net.problemParameters.m  = net.problemParameters.n - net.problemParameters.m1;
            net.problemParameters.networkSize = networkSize; %[net.problemParameters.n + net.problemParameters.m1, net.problemParameters.n + net.problemParameters.m1];
            
            net.problemParameters.P = zeros(net.problemParameters.networkSize);
            net.problemParameters.q = zeros(net.problemParameters.networkSize(1),1);
            net.problemParameters.P(1:end-net.problemParameters.m1,1:end-net.problemParameters.m1) = net.originalParameters.P;
            net.problemParameters.q(1:end-net.problemParameters.m1) = net.originalParameters.q;

            if net.problemParameters.m1 > 0
                net.problemParameters.R = zeros(net.problemParameters.m,net.problemParameters.n + net.problemParameters.m1);
                net.problemParameters.R(1:net.problemParameters.m1,1:net.problemParameters.n) = net.originalParameters.A;
                net.problemParameters.R(net.problemParameters.m1+1:net.problemParameters.m,1:net.problemParameters.n) = net.originalParameters.Aeq;
                net.problemParameters.R(1:net.problemParameters.m1,net.problemParameters.n+1:net.problemParameters.n+net.problemParameters.m1) = eye(net.problemParameters.m1);
                net.problemParameters.b = zeros(net.problemParameters.m,1);
                net.problemParameters.b(1:net.problemParameters.m1) = net.originalParameters.b;
                net.problemParameters.b(net.problemParameters.m1+1:net.problemParameters.m) = net.originalParameters.beq;
            else
                net.problemParameters.R = net.originalParameters.Aeq;
                net.problemParameters.b = net.originalParameters.beq;
            end
                        
            % Setting default options if not provided
            if nargin < 5            
                net = addDefaultOptionValues(net, options);           
                
            % Setting desired options available in options structure
            else    
                net = setOptions(net, options);
                net = addDefaultOptionValues(net, options);
            end           
            
            % Initialize results;
            net = init(net);            
        end
    end
    
    % --- Methods definitions --- %
	methods (Hidden = true, Access = private)
        net = addDefaultOptionValues(net, options);
        net = setOptions(net,options);        
    end
    
	methods (Static = true, Hidden = true, Access = protected)
        isValidSetting(property,value);        
    end
    
    methods (Static = true, Hidden = false, Access = public)
        options = createOptions(varargin);
    end
    
	methods (Hidden = false, Access = protected)
        net = init(net);
    end
end
