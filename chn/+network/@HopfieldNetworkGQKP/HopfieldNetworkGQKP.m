classdef HopfieldNetworkGQKP < network.HopfieldNetwork
    % HopfieldNetworkGQKP Generic Continuous Hopfield Network (CHN)
    %
    %   This class is used to create Generic Hopfield Networks, that is,
    %   Hopfield Networks to solve the following Generic Quadratic Knapsack 
    %   Problems (GQKP): 
    %
    %      min{1/2 * v' * P * v + q' * v}
    %      subject to
    %        A * v <= b
    %        Aeq * v = beq
    %
    %   which internally gets transformed to:
    % 
    %      min{1/2 * v' * P * v + q' * v}
    %      subject to
    %        R * v = b
    %        v_i belongs to {0,1} for all i = 1,...n
    %        v_{n+k} belongs to [0,1], k = 1,..., m1
    %
    %   HopfieldNetworkGQKP input arguments:
    %       P                          - Square matrix to define the
    %                                    minimization problem.
    %       q                          - Column vector to define the
    %                                    minimization problem.
    %       A                          - Matrix to define the 
    %                                    linear inequality constraints.
    %       b                          - Vector to define the 
    %                                    linear inequality constraints.
    %       Aeq                        - Matrix to define the 
    %                                    linear equality constraints.
    %       beq                        - Vector to define the 
    %                                    linear equality constraints.
    %       options                    - Options for simulation of the
    %                                    Hopfield Network (object from the
    %                                    HopfieldNetworkGQKPOptions class).
    % 
    %   Example:
    %       Build a Generic Hopfield Network to solve the following GQKP: 
    %
    %          min{1/2 (4*v_{1}^2 - 2*v_{2}^2)}
    %          subject to v_{1} + v_{2} = 1
    %       
    %       P = [4,0;0,-2];
    %       q = [0;0];
    %       Aeq = [1,1];
    %       beq = 1;
    %       A = [];
    %       b = [];
    %       net = network.HopfieldNetworkGQKP(P, q, A, b, Aeq, beq);
    %
    %   See also options.HopfieldNetworkOptions, hopfieldnet,
    %   hopfieldnetOptions.
    
    properties (GetAccess = public, SetAccess = public)
        Name = 'HopfieldNetworkGQKP';
        OriginalParameters;
        ProblemParameters;
    end
    
    methods
        function net = HopfieldNetworkGQKP(P, q, A, b, Aeq, beq, opts)

            if nargin < 6 || nargin > 7
                error('HopfieldNetworkGQKP:IncorrectInputArguments','Please provide proper input arguments to HopfieldNetworkGQKP: P, q, A, b, Aeq, beq (and options)');
            else                
                % asserting inputs
                assert(isa(P,'double'), 'HopfieldNetworkGQKP:invalid_datatype', 'P must be double.');
                assert(isa(q,'double'), 'HopfieldNetworkGQKP:invalid_datatype', 'q must be double.');
                assert(isa(A,'double'), 'HopfieldNetworkGQKP:invalid_datatype', 'A must be double.');
                assert(isa(b,'double'), 'HopfieldNetworkGQKP:invalid_datatype', 'b must be double.');
                assert(isa(Aeq,'double'), 'HopfieldNetworkGQKP:invalid_datatype', 'Aeq must be double.');
                assert(isa(beq,'double'), 'HopfieldNetworkGQKP:invalid_datatype', 'beq must be double.');                
            end           
                        
            if nargin < 7
                opts = options.HopfieldNetworkGQKPOptions();
            end
            
            % Creating net object from parent class
            networkSize = size(P,1) + size(A,1);
            net = net@network.HopfieldNetwork(networkSize,opts);

            net.OriginalParameters.P = P;
            net.OriginalParameters.q = q;
            net.OriginalParameters.A = A;
            net.OriginalParameters.b = b;
            net.OriginalParameters.Aeq = Aeq;
            net.OriginalParameters.beq = beq;
                       
            net.ProblemParameters.n  = size(net.OriginalParameters.P,1);
            net.ProblemParameters.m1 = size(net.OriginalParameters.A,1);
            net.ProblemParameters.m  = net.ProblemParameters.n - net.ProblemParameters.m1;
            net.ProblemParameters.networkSize = networkSize; %[net.ProblemParameters.n + net.ProblemParameters.m1, net.ProblemParameters.n + net.ProblemParameters.m1];
            
            net.ProblemParameters.P = zeros(net.ProblemParameters.networkSize);
            net.ProblemParameters.q = zeros(net.ProblemParameters.networkSize(1),1);
            net.ProblemParameters.P(1:end-net.ProblemParameters.m1,1:end-net.ProblemParameters.m1) = net.OriginalParameters.P;
            net.ProblemParameters.q(1:end-net.ProblemParameters.m1) = net.OriginalParameters.q;

            if net.ProblemParameters.m1 > 0
                net.ProblemParameters.R = zeros(net.ProblemParameters.m,net.ProblemParameters.n + net.ProblemParameters.m1);
                net.ProblemParameters.R(1:net.ProblemParameters.m1,1:net.ProblemParameters.n) = net.OriginalParameters.A;
                net.ProblemParameters.R(net.ProblemParameters.m1+1:net.ProblemParameters.m,1:net.ProblemParameters.n) = net.OriginalParameters.Aeq;
                net.ProblemParameters.R(1:net.ProblemParameters.m1,net.ProblemParameters.n+1:net.ProblemParameters.n+net.ProblemParameters.m1) = eye(net.ProblemParameters.m1);
                net.ProblemParameters.b = zeros(net.ProblemParameters.m,1);
                net.ProblemParameters.b(1:net.ProblemParameters.m1) = net.OriginalParameters.b;
                net.ProblemParameters.b(net.ProblemParameters.m1+1:net.ProblemParameters.m) = net.OriginalParameters.beq;
            else
                net.ProblemParameters.R = net.OriginalParameters.Aeq;
                net.ProblemParameters.b = net.OriginalParameters.beq;
            end
                        
            net = setOptions(net, opts);
            
            % Initialize results;
            net = init(net);            
        end
    end
    
    % --- Methods definitions --- %
	methods (Hidden = true, Access = private)
        net = setOptions(net,options);        
    end
    
	methods (Static = true, Hidden = true, Access = private)

    end
	methods (Hidden = false, Access = protected)
        net = init(net);
    end
end
