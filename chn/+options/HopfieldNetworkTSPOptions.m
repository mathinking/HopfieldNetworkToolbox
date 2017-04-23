classdef HopfieldNetworkTSPOptions < options.HopfieldNetworkOptions
    %HopfieldNetworkTSPOptions Options for the tsphopfieldnet class. 
    %
    %   This class defines allows to define options for Hopfield Networks
    %   to solve the Traveling Salesman Problem (TSP). Then it can be used 
    %   to instantiate HopfieldNetworkTSP objects.
    %
    %   HopfieldNetworkTSPOptions properties:
    %       Cities                     - Defition of the TSP problem.
    %       Setting                    - Settings used in the Simulation of
    %                                    the Hopfield.
    %       SimFcn                     - Simulation Function.
    %       TrainFcn                   - Training Function.
    %       TrainingParam              - Training Parameters.
    %
    %   Example:
    %       Create a set of options for tsphopfieldnet to solves a TSP
    %       with 10 cities, where city coordinates are chosen randomly,
    %       distance is euclidean, maximum number of iterations is 1000, 
    %       the execution environment is the cpu and the transfer function 
    %       is a saturating linear transfer function.
    %
    %       options = tsphopfieldnetOptions('Coordinates', rand(10,2), ...
    %                                       'DistanceType', 'euc', ...
    %                                       'MaxIter', 1000, ...
    %                                       'ExecutionEnvironment', 'cpu', ...
    %                                       'TransferFcn', 'satlin');
    % 
    %   See also tsphopfieldnet, options.HopfieldNetworkOptions.
    
    properties (SetAccess = private)
        % Cities   Defition of the TSP problem.
        %   The Cities options allow to define the TSP problem through the
        %   following properties: 
        %      - Coordinates
        %      - DistanceMatrix
        %      - DistanceType
        %      - Names
        %      - Subtours
        %      - SubtoursPositions
        %      - Tau
        Cities
    end
    properties (SetAccess = private) 
        % SimFcn   Simulation Function.
        %   The Simulation Function determines how simulation takes place
        %   in the Hopfield Network. Possible values: 
        %
        %      'euler', 'talavan-yanez' (default) and 'divide-conquer'
        SimFcn
        
        % TrainFcn   Training Function.
        %   The Training Function determines how training takes place in
        %   the Hopfield Network. Possible values:
        %
        %      'trainty' (default)
        TrainFcn
        
        % TrainParam   Training Parameters.
        %   The Training Parameters are used to determine the appropriate
        %   weigths and biases of the Hopfield Network. 
        TrainParam
    end
    
    methods
        function opts = HopfieldNetworkTSPOptions(varargin)
            % Parsing options to be set by hopfieldnetworkOptions class
            opts = opts@options.HopfieldNetworkOptions(varargin{:});
            
            parser = inputParser;
            
            defaultTrainFcn = 'trainty';
            defaultSimFcn   = 'talavan-yanez'; 
            
            defaultTrainParamK = 0;
            
            defaultCitiesNames = '';
            defaultCitiesDistanceType = 'EUC';
            defaultCitiesCoordinates = []; 
            defaultCitiesDistanceMatrix = [];
            defaultCitiesSubtours = '';
            defaultCitieSubtoursPositions = [];
            defaultCitiesTau = [];
            
            parser.addParameter('TrainFcn', defaultTrainFcn, @opts.iIsChar);
            parser.addParameter('SimFcn', defaultSimFcn, @opts.iIsChar);
            parser.addParameter('K', defaultTrainParamK, @opts.iIsIntegerScalarGreaterOrEqualThanZero);
            parser.addParameter('Names', defaultCitiesNames, @opts.iIsCellStr)
            parser.addParameter('DistanceType', defaultCitiesDistanceType, @opts.iIsChar)
            parser.addParameter('Coordinates', defaultCitiesCoordinates, @opts.iIsRealNumericWithTwoColumns)
            parser.addParameter('DistanceMatrix', defaultCitiesDistanceMatrix, @opts.iIsRealNumericSquareMatrix)
            parser.addParameter('Subtours', defaultCitiesSubtours, @opts.iIsCellStr)
            parser.addParameter('SubtoursPositions', defaultCitieSubtoursPositions, @opts.iIsInteger)
            parser.addParameter('Tau', defaultCitiesTau, @opts.iIsIntegerScalarGreaterOrEqualThanZero)
            
            % Parsing remaining options
            parser.parse(opts.Unmatched);
            
            % Clearing Unmatched from memory - any Unmatched property
            % should have produced an error
            opts.Unmatched = [];
            
            opts.TrainFcn = opts.iMatchWithValidTrainFcn(parser.Results.TrainFcn);
            opts.SimFcn   = opts.iMatchWithValidSimFcn(parser.Results.SimFcn);
            opts.TrainParam.K = parser.Results.K;
            opts.Cities.Names = parser.Results.Names;
            opts.Cities.DistanceType = opts.iMatchWithValidDistanceType(parser.Results.DistanceType);
            opts.Cities.Coordinates = parser.Results.Coordinates;
            opts.Cities.DistanceMatrix = opts.iVerifyExistenceIfExplicit(parser.Results.DistanceMatrix, opts.Cities.DistanceType, opts.Cities.Coordinates);
            opts.Cities.Subtours = parser.Results.Subtours;
            opts.Cities.SubtoursPositions = opts.iMatchWithSubtoursLength(parser.Results.SubtoursPositions, opts.Cities.Subtours);
            opts.Cities.Tau = parser.Results.Tau;

            opts.Cities = orderfields(opts.Cities);            
        end
    end
    methods (Static = true, Access = private)
        function tf = iIsChar(x)
            tf = ischar(x);
        end
        function tf = iIsIntegerScalarGreaterOrEqualThanZero(x)
            tf = isscalar(x) && options.HopfieldNetworkTSPOptions.iIsInteger(x) && (x >= 0);
        end
        function tf = iIsInteger(x)
            tf = options.HopfieldNetworkTSPOptions.iIsRealNumeric(x) && all(mod(x,1)==0);
        end
        function tf = iIsRealNumeric(x)
            tf = isreal(x) && isnumeric(x);
        end
        function tf = iIsCellStr(x)
            tf = iscell(x) && all(cellfun(@ischar,x));
        end
        function tf = iIsRealNumericWithTwoColumns(x)
            tf = options.HopfieldNetworkTSPOptions.iIsRealNumeric(x) && (size(x,2) == 2 || isempty(x));
        end
        function tf = iIsRealNumericSquareMatrix(x)
            tf = options.HopfieldNetworkTSPOptions.iIsRealNumeric(x) && (size(x,1) == size(x,2));
        end

        function chosenTrainFcn = iMatchWithValidTrainFcn(inputTrainFcn)
            validTrainFcn = {'trainty'};
            chosenTrainFcn = validatestring(inputTrainFcn, validTrainFcn, 'HopfieldNetworkTSPOptions', 'TrainFcn');
        end
        function chosenSimFcn = iMatchWithValidSimFcn(inputSimFcn)
            validSimFcn = {'euler','talavan-yanez','divide-conquer'};
            chosenSimFcn = validatestring(inputSimFcn, validSimFcn, 'HopfieldNetworkTSPOptions', 'SimFcn');
        end         
        function chosenType = iMatchWithValidDistanceType(inputType)
            validType = {'geo','euc_2d','euc','att','ceil_2d','explicit'};
            chosenType = validatestring(inputType, validType, 'HopfieldNetworkTSPOptions', 'DistanceType');
        end    
        function inputSubtoursPositions = iMatchWithSubtoursLength(inputSubtoursPositions, inputSubtours)
            assert(length(inputSubtoursPositions) == length(inputSubtours), 'HopfieldNetworkTSPOptions:unevensubtours', ...
                'The number of Subtours starting positions does not match the number of Subtours.');
        end
        function inputDistanceMatrix = iVerifyExistenceIfExplicit(inputDistanceMatrix, inputDistanceType, inputCoordinates)
            if ~isempty(inputDistanceMatrix) && isempty(inputCoordinates)
                assert(strcmp(inputDistanceType, 'explicit'), 'HopfieldNetworkTSPOptions:DistanceMatrixGivenTypeNotExplicit', ...
                    'Set ''DistanceType'' to ''explicit'' if ''DistanceMatrix'' is provided and ''Coordinates'' are empty.');
            end
            if strcmp(inputDistanceType, 'explicit')
                assert(~isempty(inputDistanceMatrix), 'HopfieldNetworkTSPOptions:emptyDistanceMatrix', ...
                    '''DistanceMatrix'' cannot be empty if ''DistanceType'' is set to ''explicit''.');
                assert(isempty(inputCoordinates), 'HopfieldNetworkTSPOptions:nonEmptyCoordinates', ...
                    '''Coordinates'' must be empty if ''DistanceType'' is set to ''explicit''.'); 
            end
        end       
    end
    
end
