classdef HopfieldNetworkOptions < handle
    %HopfieldNetworkOptions Options for the HopfieldNetwork class. 
    %
    %   This class defines abstract methods and properties to set options 
    %   for Hopfield Networks. Therefore it can only be instantiated from 
    %   within a child class.
    %
    %   HopfieldNetworkOptions properties:
    %       Setting                    - Settings used in the Simulation of
    %                                    the Hopfield Network.
    %       SimFcn                     - Simulation Function.
    %       TrainFcn                   - Training Function.
    %       TrainingParam              - Training Parameters.
    %       Unmatched                  - Unmatched Properties.
    %
    %   See also options.HopfieldNetworkGQKPOptions,
    %   options.HopfieldNetworkTSPOptions.
    
    properties (SetAccess = private)
        % Setting   Settings used in the Simulation of the Hopfield Network. 
        %   The following Settings can be set for this class, though they 
        %   can only be used in a child class, hopfieldnet or 
        %   tsphopfieldnet:
        %      - CheckpointPath
        %      - Dt
        %      - E
        %      - ExecutionEnvironment
        %      - MaxIter
        %      - Q
        %      - R_Iter
        %      - SimulationPlot
        %      - SimulationPlotPauseTime
        %      - TransferFcn
        %      - U0
        %      - Verbose
        %      - VerboseFrequency
        Setting
    end
    properties (SetAccess = private, Abstract = true)
        % Scheme   Network scheme
        %   The Scheme determines the network pipeline used in the
        %   simulation of the Continuous Hopfield Network
        Scheme
        
        % SimFcn   Simulation Function.
        %   The Simulation Function determines how simulation takes place
        %   in the Hopfield Network. This property is used in the sim 
        %   function of its deriving class: hopfieldnetOptions or 
        %   tsphopfieldnetOptions.
        SimFcn
        
        % TrainFcn   Training Function.
        %   The Training Function determines how training takes place in
        %   the Hopfield Network and therefore how weights and biases of
        %   the network are set. This property is implemented in its
        %   deriving class: hopfieldnetOptions or tsphopfieldnetOptions.
        TrainFcn

        % TrainParam   Training Parameters.
        %   The Training Parameters are used to determine the appropriate
        %   weigths and biases of the Hopfield Network. This property is 
        %   implemented in its deriving class: hopfieldnetOptions or 
        %   tsphopfieldnetOptions.       
        TrainParam
    end
    properties (GetAccess = protected, SetAccess = protected)
        % Unmatched   Unmatched Properties.
        %   Unmatched properties may be set in the hopfieldnetworkOptions
        %   child classes: hopfieldnetOptions or tsphopfieldnetOptions.       
        Unmatched
    end
    
    methods
        function opts = HopfieldNetworkOptions(varargin)
            parser = inputParser;
            parser.KeepUnmatched = true;
            
            defaultSettingU0 = 0.3;
            defaultSettingTransferFcn = 'tanh';
            defaultSettingVerbose = false;
            defaultSettingVerboseFrequency = NaN;
            defaultSettingExecutionEnvironment = 'CPU';
            defaultSettingMaxIter = 2000;
            defaultSettingE = 13;
            defaultSettingQ = 0.8;
            defaultSettingR_Iter = 20;
            defaultSettingDt = 0.01;
            defaultSettingCheckpointPath = '';
            defaultSettingSimulationPlot = false;
            defaultSettingSimulationPlotPauseTime = 0.3; 

            parser.addParameter('U0', defaultSettingU0, @opts.iIsRealNumericScalarGreaterThanZero);
            parser.addParameter('TransferFcn', defaultSettingTransferFcn);
            parser.addParameter('Verbose', defaultSettingVerbose, @opts.iIsScalarAndLogicalOneOrZero);
            parser.addParameter('VerboseFrequency', defaultSettingVerboseFrequency, @opts.iIsPositiveIntegerScalar);
            parser.addParameter('ExecutionEnvironment', defaultSettingExecutionEnvironment);
            parser.addParameter('MaxIter', defaultSettingMaxIter, @opts.iIsPositiveIntegerScalar);
            parser.addParameter('E', defaultSettingE, @opts.iIsPositiveIntegerScalar);
            parser.addParameter('Q', defaultSettingQ, @opts.iIsRealNumericScalarGreaterThanZeroAndLessOrEqualThanOne);
            parser.addParameter('R_Iter', defaultSettingR_Iter, @opts.iIsPositiveIntegerScalar);
            parser.addParameter('Dt', defaultSettingDt, @opts.iIsRealNumericScalarGreaterThanZero);
            parser.addParameter('CheckpointPath', defaultSettingCheckpointPath, @opts.iIsValidCheckpointPath);
            parser.addParameter('SimulationPlot', defaultSettingSimulationPlot, @opts.iIsScalarAndLogicalOneOrZero);
            parser.addParameter('SimulationPlotPauseTime', defaultSettingSimulationPlotPauseTime, @opts.iIsRealNumericScalarGreaterThanZero);

            parser.parse(varargin{:});

            opts.Setting.U0 = parser.Results.U0;
            opts.Setting.TransferFcn = opts.iMatchWithValidTransferFunction(parser.Results.TransferFcn);
            opts.Setting.Verbose = logical(parser.Results.Verbose);
            opts.Setting.VerboseFrequency = parser.Results.VerboseFrequency;
            opts.Setting.ExecutionEnvironment = opts.iMatchWithValidExecutionEnvironment(parser.Results.ExecutionEnvironment);
            opts.Setting.MaxIter = parser.Results.MaxIter;
            opts.Setting.E = parser.Results.E;
            opts.Setting.Q = parser.Results.Q;
            opts.Setting.R_Iter = parser.Results.R_Iter;
            opts.Setting.Dt = parser.Results.Dt;
            opts.Setting.CheckpointPath = parser.Results.CheckpointPath;
            opts.Setting.SimulationPlot = parser.Results.SimulationPlot;
            opts.Setting.SimulationPlotPauseTime = parser.Results.SimulationPlotPauseTime;
            
            opts.Setting = orderfields(opts.Setting);
            opts.Unmatched = parser.Unmatched;
        end
    end
    methods (Static = true, Access = private)
        function tf = iIsRealNumericScalar(x)
            tf = isscalar(x) && isreal(x) && isnumeric(x);
        end
        function tf = iIsRealNumericScalarGreaterThanZero(x)
            tf = options.HopfieldNetworkOptions.iIsRealNumericScalar(x) && (x > 0);
        end
        function tf = iIsScalarAndLogicalOneOrZero(x)
            tf = isscalar(x) && options.HopfieldNetworkOptions.iIsLogicalOneOrZero(x);
        end
        function tf =  iIsLogicalOneOrZero(x)
            tf = islogical(x) || (x == 1) || (x == 0);
        end
        function tf = iIsPositiveIntegerScalar(x)
            tf = isscalar(x) && options.HopfieldNetworkOptions.iIsInteger(x) && (x > 0);
        end
        function tf = iIsInteger(x)
            tf = isreal(x) && isnumeric(x) && all(mod(x,1)==0);
        end
        function tf = iIsRealNumericScalarGreaterThanZeroAndLessOrEqualThanOne(x)
            tf = options.HopfieldNetworkOptions.iIsRealNumericScalar(x) && (x > 0) && (x <= 1);
        end
        function tf = iIsValidCheckpointPath(x)
        % iIsValidCheckpointPath Return true if x is a valid checkpoint path or
        % an empty string. Valid checkpoint paths are existing directories with
        % write access.
            if ~ischar(x)
                error('HopfieldNetworkOptions:notChar', '''CheckpointPath'' must be a character array.');
            elseif isempty(x)
                tf = true;
            elseif ~isdir(x) 
                error('HopfieldNetworkOptions:notExistingPath', 'Provide an existing path for ''CheckpointPath''.');
            elseif ~options.HopfieldNetworkOptions.iCanWriteToDir(x)
                error('HopfieldNetworkOptions:notWritablePath', 'Provide a ''CheckpointPath'' with write access.');
            else
                tf = true;
            end
        end
        function tf = iCanWriteToDir(proposedDir)
            [~, status] = fileattrib(proposedDir);
            tf = status.UserWrite;
        end
        function chosenExecutionEnvironment = iMatchWithValidExecutionEnvironment(inputEnvironment)
            validExecutionEnvironments = {'cpu', 'gpu'};
            chosenExecutionEnvironment = validatestring(inputEnvironment, validExecutionEnvironments, 'HopfieldNetworkOptions', 'ExecutionEnvironment');
        end
        function chosenTransferFnc = iMatchWithValidTransferFunction(inputTransferFcn)
            validTransferFcns = {'satlin', 'tanh'};
            chosenTransferFnc = validatestring(inputTransferFcn, validTransferFcns, 'HopfieldNetworkOptions', 'TransferFcn');
        end
    end
    
end
