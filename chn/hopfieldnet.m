function net = hopfieldnet(P, q, A, b, Aeq, beq, options)
% hopfieldnet   Generic Continuous Hopfield Network (CHN)
%
%   net = hopfieldnet(P, q, A, b, Aeq, beq) creates a Hopfield
%   Network to solve the Generic Quadratic Knapsack 
%   Problems (GQKP) using a default set of options.
% 
%   net = hopfieldnet(P, q, A, b, Aeq, beq, options) specifies optional
%   parameters that allow to configure the Hopfield Network.
%
%   A Generic Quadratic Knapsack Problem (GQKP) is defined in the
%   following way
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
%   hopfieldnet input arguments:
%
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
%       options                    - (optional) Options for simulation 
%                                    of the Hopfield Network (object 
%                                    from the HopfieldNetworkGQKPOptions
%                                    class). This options object can be
%                                    created using the
%                                    hopfieldnetOptions function.
%
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
%       net = hopfieldnet(P, q, A, b, Aeq, beq);
%
%   See also hopfieldnetOptions, options.HopfieldNetworkOptions.

    if nargin == 6
        net = network.HopfieldNetworkGQKP(P, q, A, b, Aeq, beq);
    elseif nargin == 7
        net = network.HopfieldNetworkGQKP(P, q, A, b, Aeq, beq, options);
    else
        error('HopfieldNetworkGQKP:IncorrectInputArguments','Please provide proper input arguments to HopfieldNetworkGQKP: P, q, A, b, Aeq, beq (and options)');
    end
    
end
