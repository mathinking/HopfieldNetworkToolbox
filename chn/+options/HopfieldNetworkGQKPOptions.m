classdef HopfieldNetworkGQKPOptions < options.HopfieldNetworkOptions
    %HopfieldNetworkGQKPOptions Options for the HopfieldNetworkGQKP class. 
    %
    %   This class defines allows to define options for Generic Hopfield 
    %   Networks. Then it can be used to instantiate HopfieldNetworkGQKP
    %   objects.
    %
    %   HopfieldNetworkGQKPOptions properties:
    %       Setting                    - Settings used in the Simulation of
    %                                    the Hopfield.
    %       SimFcn                     - Simulation Function.
    %       TrainFcn                   - Training Function.
    %       TrainingParam              - Training Parameters.
    %
    %   Example:
    %       Create a set of options for HopfieldNetworkGQKP to solve a 
    %       Generic Hopfield Network Problem, where the simulation function 
    %       is using euler's method, maximum number of iterations is 5000 
    %       and the transfer function is a hyperbolic tangent transfer
    %       function.
    %       
    %       options = options.HopfieldNetworkGQKPOptions('SimFcn', 'euler', ...
    %                                                    'MaxIter', 5000, ...
    %                                                    'TransferFcn', 'tanh');
    % 
    %   See also hopfieldnet, options.HopfieldNetworkOptions.
    
    properties (SetAccess = private)
        % Scheme   Network scheme
        %   The Scheme determines the network pipeline used in the
        %   simulation of the Continuous Hopfield Network
        Scheme = 'classic'
        
        % SimFcn   Simulation Function.
        %   The Simulation Function determines how simulation takes place
        %   in the Hopfield Network. Possible values: 
        %
        %      'euler' and 'talavan-yanez' (default)
        SimFcn
        
        % TrainFcn   Training Function.
        %   The Training Function determines how training takes place in
        %   the Hopfield Network. Possible values:
        %
        %      'traingty' (default)
        TrainFcn
        
        % TrainParam   Training Parameters.
        %   The Training Parameters are used to determine the appropriate
        %   weigths and biases of the Hopfield Network. 
        TrainParam
    end
    
    methods
        function opts = HopfieldNetworkGQKPOptions(varargin)
            opts = opts@options.HopfieldNetworkOptions(varargin{:});
            
            parser = inputParser;

            defaultTrainFcn = 'traingty';
            defaultSimFcn   = 'talavan-yanez'; 

            parser.addParameter('TrainFcn', defaultTrainFcn, @opts.iIsChar);
            parser.addParameter('SimFcn', defaultSimFcn, @opts.iIsChar);
            
            parser.parse(opts.Unmatched);
            
            % Clearing Unmatched from memory - any Unmatched property
            % should have produced an error
            opts.Unmatched = [];
            
            opts.SimFcn   = opts.iMatchWithValidSimFcn(parser.Results.SimFcn);
            opts.TrainFcn = opts.iMatchWithValidTrainFcn(parser.Results.TrainFcn); 
        end
    end    
    methods (Static = true, Access = private)
        function tf = iIsChar(x)
            tf = ischar(x);
        end
        function chosenTrainFcn = iMatchWithValidTrainFcn(inputTrainFcn)
            validTrainFcn = {'traingty'};
            chosenTrainFcn = validatestring(inputTrainFcn, validTrainFcn, 'HopfieldNetworkGQKPOptions', 'TrainFcn');
        end
        function chosenSimFcn = iMatchWithValidSimFcn(inputSimFcn)
            validSimFcn = {'euler','runge-kutta','talavan-yanez'};
            chosenSimFcn = validatestring(inputSimFcn, validSimFcn, 'HopfieldNetworkGQKPOptions', 'SimFcn');
        end          
    end

end
