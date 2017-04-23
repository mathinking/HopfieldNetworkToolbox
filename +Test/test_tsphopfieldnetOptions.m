classdef test_tsphopfieldnetOptions < matlab.unittest.TestCase
%TEST_TSPHOPFIELDNETOPTIONS Summary of this class goes here
%   Detailed explanation goes here
    
	methods (Test)
        % [the name of the tested method]_[expected input / tested state]_[expected behavior]
        
        % Verify inputs
        % Checking tsphopfieldnetOptions input arguments
        function HopfieldNetworkOptions_U0NotGreaterThanZero_Errors(testCase)
            verifyError(testCase,@()tsphopfieldnetOptions('U0',0),'MATLAB:InputParser:ArgumentFailedValidation');
        end
        function HopfieldNetworkOptions_U0NotScalar_Errors(testCase)
            verifyError(testCase,@()tsphopfieldnetOptions('U0',[1,1]),'MATLAB:InputParser:ArgumentFailedValidation');
        end
        function HopfieldNetworkOptions_U0NotReal_Errors(testCase)
            verifyError(testCase,@()tsphopfieldnetOptions('U0',1+2i),'MATLAB:InputParser:ArgumentFailedValidation');
        end
        function HopfieldNetworkOptions_U0NotNumeric_Errors(testCase)
            verifyError(testCase,@()tsphopfieldnetOptions('U0','1'),'MATLAB:InputParser:ArgumentFailedValidation');
        end
        function HopfieldNetworkOptions_TransferFcnNotValid_Errors(testCase)
            verifyError(testCase,@()tsphopfieldnetOptions('TransferFcn','linear'),'MATLAB:HopfieldNetworkOptions:unrecognizedStringChoice');
        end
        function HopfieldNetworkOptions_VerboseNotScalar_Errors(testCase)
            verifyError(testCase,@()tsphopfieldnetOptions('Verbose',[true, false]),'MATLAB:InputParser:ArgumentFailedValidation');
        end
        function HopfieldNetworkOptions_VerboseNotLogical_Errors(testCase)
            verifyError(testCase,@()tsphopfieldnetOptions('Verbose',5),'MATLAB:InputParser:ArgumentFailedValidation');
        end
        function HopfieldNetworkOptions_VerboseFrequencyNotScalar_Errors(testCase)
            verifyError(testCase,@()tsphopfieldnetOptions('VerboseFrequency',[10,10]),'MATLAB:InputParser:ArgumentFailedValidation');
        end
        function HopfieldNetworkOptions_VerboseFrequencyNotReal_Errors(testCase)
            verifyError(testCase,@()tsphopfieldnetOptions('VerboseFrequency',1+2i),'MATLAB:InputParser:ArgumentFailedValidation');
        end
        function HopfieldNetworkOptions_VerboseFrequencyNotNumeric_Errors(testCase)
            verifyError(testCase,@()tsphopfieldnetOptions('VerboseFrequency','freq'),'MATLAB:InputParser:ArgumentFailedValidation');
        end
        function HopfieldNetworkOptions_VerboseFrequencyNotInteger_Errors(testCase)
            verifyError(testCase,@()tsphopfieldnetOptions('VerboseFrequency',50.5),'MATLAB:InputParser:ArgumentFailedValidation');
        end
        function HopfieldNetworkOptions_VerboseFrequencyNotGreaterZero_Errors(testCase)
            verifyError(testCase,@()tsphopfieldnetOptions('VerboseFrequency',0),'MATLAB:InputParser:ArgumentFailedValidation');
        end
        function HopfieldNetworkOptions_ExecutionEnvironmentNotValid_Errors(testCase)
            verifyError(testCase,@()tsphopfieldnetOptions('ExecutionEnvironment','fpga'),'MATLAB:HopfieldNetworkOptions:unrecognizedStringChoice');
        end
        function HopfieldNetworkOptions_MaxIterNotScalar_Errors(testCase)
            verifyError(testCase,@()tsphopfieldnetOptions('MaxIter',[2000,2000]),'MATLAB:InputParser:ArgumentFailedValidation');           
        end
        function HopfieldNetworkOptions_MaxIterNotReal_Errors(testCase)
            verifyError(testCase,@()tsphopfieldnetOptions('MaxIter',1+2i),'MATLAB:InputParser:ArgumentFailedValidation');           
        end
        function HopfieldNetworkOptions_MaxIterNotNumeric_Errors(testCase)
            verifyError(testCase,@()tsphopfieldnetOptions('MaxIter','iter'),'MATLAB:InputParser:ArgumentFailedValidation');           
        end
        function HopfieldNetworkOptions_MaxIterNotInteger_Errors(testCase)
            verifyError(testCase,@()tsphopfieldnetOptions('MaxIter',2000.1),'MATLAB:InputParser:ArgumentFailedValidation');           
        end
        function HopfieldNetworkOptions_MaxIterNotGreaterThanZero_Errors(testCase)
            verifyError(testCase,@()tsphopfieldnetOptions('MaxIter',0),'MATLAB:InputParser:ArgumentFailedValidation');           
        end
        function HopfieldNetworkOptions_ENotScalar_Errors(testCase)
            verifyError(testCase,@()tsphopfieldnetOptions('E',[13,13]),'MATLAB:InputParser:ArgumentFailedValidation');           
        end
        function HopfieldNetworkOptions_ENotReal_Errors(testCase)
            verifyError(testCase,@()tsphopfieldnetOptions('E',1+2i),'MATLAB:InputParser:ArgumentFailedValidation');           
        end
        function HopfieldNetworkOptions_ENotNumeric_Errors(testCase)
            verifyError(testCase,@()tsphopfieldnetOptions('E','e'),'MATLAB:InputParser:ArgumentFailedValidation');           
        end
        function HopfieldNetworkOptions_ENotInteger_Errors(testCase)
            verifyError(testCase,@()tsphopfieldnetOptions('E',12.1),'MATLAB:InputParser:ArgumentFailedValidation');           
        end
        function HopfieldNetworkOptions_ENotGreaterThanZero_Errors(testCase)
            verifyError(testCase,@()tsphopfieldnetOptions('E',0),'MATLAB:InputParser:ArgumentFailedValidation');           
        end
        function HopfieldNetworkOptions_QNotScalar_Errors(testCase)
            verifyError(testCase,@()tsphopfieldnetOptions('Q',[13,13]),'MATLAB:InputParser:ArgumentFailedValidation');           
        end
        function HopfieldNetworkOptions_QNotReal_Errors(testCase)
            verifyError(testCase,@()tsphopfieldnetOptions('Q',1+2i),'MATLAB:InputParser:ArgumentFailedValidation');           
        end
        function HopfieldNetworkOptions_QNotNumeric_Errors(testCase)
            verifyError(testCase,@()tsphopfieldnetOptions('Q','e'),'MATLAB:InputParser:ArgumentFailedValidation');           
        end
        function HopfieldNetworkOptions_QNotGreaterThanZero_Errors(testCase)
            verifyError(testCase,@()tsphopfieldnetOptions('Q',0),'MATLAB:InputParser:ArgumentFailedValidation');           
        end
        function HopfieldNetworkOptions_QNotLessOrEqualThanOne_Errors(testCase)
            verifyError(testCase,@()tsphopfieldnetOptions('Q',1+eps),'MATLAB:InputParser:ArgumentFailedValidation');           
        end        
        function HopfieldNetworkOptions_RIterNotScalar_Errors(testCase)
            verifyError(testCase,@()tsphopfieldnetOptions('R_Iter',[20,20]),'MATLAB:InputParser:ArgumentFailedValidation');           
        end
        function HopfieldNetworkOptions_RIterNotReal_Errors(testCase)
            verifyError(testCase,@()tsphopfieldnetOptions('R_Iter',1+2i),'MATLAB:InputParser:ArgumentFailedValidation');           
        end
        function HopfieldNetworkOptions_RIterNotNumeric_Errors(testCase)
            verifyError(testCase,@()tsphopfieldnetOptions('R_Iter','R_ITER'),'MATLAB:InputParser:ArgumentFailedValidation');           
        end
        function HopfieldNetworkOptions_RIterNotInteger_Errors(testCase)
            verifyError(testCase,@()tsphopfieldnetOptions('R_Iter',20.1),'MATLAB:InputParser:ArgumentFailedValidation');           
        end
        function HopfieldNetworkOptions_RIterNotGreaterThanZero_Errors(testCase)
            verifyError(testCase,@()tsphopfieldnetOptions('R_Iter',0),'MATLAB:InputParser:ArgumentFailedValidation');           
        end
        function HopfieldNetworkOptions_DtNotScalar_Errors(testCase)
            verifyError(testCase,@()tsphopfieldnetOptions('Dt',[0.001,0.01]),'MATLAB:InputParser:ArgumentFailedValidation');           
        end
        function HopfieldNetworkOptions_DtNotReal_Errors(testCase)
            verifyError(testCase,@()tsphopfieldnetOptions('Dt',1i),'MATLAB:InputParser:ArgumentFailedValidation');           
        end
        function HopfieldNetworkOptions_DtNotNumeric_Errors(testCase)
            verifyError(testCase,@()tsphopfieldnetOptions('Dt','small'),'MATLAB:InputParser:ArgumentFailedValidation');           
        end
        function HopfieldNetworkOptions_DtNotGreaterThanZero_Errors(testCase)
            verifyError(testCase,@()tsphopfieldnetOptions('Dt',0'),'MATLAB:InputParser:ArgumentFailedValidation');           
        end
        function HopfieldNetworkOptions_CheckpointPathNotChar_Errors(testCase)
            verifyError(testCase,@()tsphopfieldnetOptions('CheckpointPath',0),'HopfieldNetworkOptions:notChar');
        end
        function HopfieldNetworkOptions_CheckpointPathNotDirectory_Errors(testCase)
            verifyError(testCase,@()tsphopfieldnetOptions('CheckpointPath','C:/madeUpDirectory'),'HopfieldNetworkOptions:notExistingPath');
        end
        function HopfieldNetworkOptions_CheckpointPathNotWritable_Errors(testCase)
            nonWritablePath = fullfile(pwd,'CheckpointPathNotWritable');
            mkdir(nonWritablePath);
            fileattrib(nonWritablePath,'-w')
            verifyError(testCase,@()tsphopfieldnetOptions('CheckpointPath',nonWritablePath),'HopfieldNetworkOptions:notWritablePath');
            fileattrib(nonWritablePath,'+w')
            rmdir(nonWritablePath);
        end
        function HopfieldNetworkOptions_SimulationPlotNotScalar_Errors(testCase)
            verifyError(testCase,@()tsphopfieldnetOptions('SimulationPlot',[true, false]),'MATLAB:InputParser:ArgumentFailedValidation');
        end
        function HopfieldNetworkOptions_SimulationPlotNotLogical_Errors(testCase)
            verifyError(testCase,@()tsphopfieldnetOptions('SimulationPlot',5),'MATLAB:InputParser:ArgumentFailedValidation');
        end
        function HopfieldNetworkOptions_SimPlotPauseTimeNotGreaterZero_Errors(testCase)
            verifyError(testCase,@()tsphopfieldnetOptions('SimulationPlotPauseTime',0),'MATLAB:InputParser:ArgumentFailedValidation');
        end
        function HopfieldNetworkOptions_SimPlotPauseTimeNotScalar_Errors(testCase)
            verifyError(testCase,@()tsphopfieldnetOptions('SimulationPlotPauseTime',[1,1]),'MATLAB:InputParser:ArgumentFailedValidation');
        end
        function HopfieldNetworkOptions_SimPlotPauseTimeNotReal_Errors(testCase)
            verifyError(testCase,@()tsphopfieldnetOptions('SimulationPlotPauseTime',1+2i),'MATLAB:InputParser:ArgumentFailedValidation');
        end
        function HopfieldNetworkOptions_SimPlotPauseTimeNotNumeric_Errors(testCase)
            verifyError(testCase,@()tsphopfieldnetOptions('SimulationPlotPauseTime','1'),'MATLAB:InputParser:ArgumentFailedValidation');
        end
        function tsphopfieldnetOptions_TrainFcnNotChar_Errors(testCase)
            verifyError(testCase,@()tsphopfieldnetOptions('TrainFcn',0),'MATLAB:InputParser:ArgumentFailedValidation');
        end        
        function tsphopfieldnetOptions_TrainFcnNotMatchValidOptions_Errors(testCase)
            verifyError(testCase,@()tsphopfieldnetOptions('TrainFcn','traingty'),'MATLAB:HopfieldNetworkTSPOptions:unrecognizedStringChoice');
        end
        function tsphopfieldnetOptions_SimFcnNotChar_Errors(testCase)
            verifyError(testCase,@()tsphopfieldnetOptions('SimFcn',0),'MATLAB:InputParser:ArgumentFailedValidation');
        end        
        function tsphopfieldnetOptions_SimFcnNotMatchValidOptions_Errors(testCase)
            verifyError(testCase,@()tsphopfieldnetOptions('SimFcn','yanez-talavan'),'MATLAB:HopfieldNetworkTSPOptions:unrecognizedStringChoice');
        end
        function tsphopfieldnetOptions_KNotScalar_Errors(testCase)
            verifyError(testCase,@()tsphopfieldnetOptions('K',[10,10]),'MATLAB:InputParser:ArgumentFailedValidation');
        end
        function tsphopfieldnetOptions_KNotReal_Errors(testCase)
            verifyError(testCase,@()tsphopfieldnetOptions('K',1+2i),'MATLAB:InputParser:ArgumentFailedValidation');           
        end
        function tsphopfieldnetOptions_KNotNumeric_Errors(testCase)
            verifyError(testCase,@()tsphopfieldnetOptions('K','5'),'MATLAB:InputParser:ArgumentFailedValidation');           
        end       
        function tsphopfieldnetOptions_KNotInteger_Errors(testCase)
            verifyError(testCase,@()tsphopfieldnetOptions('K',10.1),'MATLAB:InputParser:ArgumentFailedValidation');           
        end
        function tsphopfieldnetOptions_KNotGreaterOrEqualThanZero_Errors(testCase)
            verifyError(testCase,@()tsphopfieldnetOptions('K',-1),'MATLAB:InputParser:ArgumentFailedValidation');           
        end
        function tsphopfieldnetOptions_NamesNotCell_Errors(testCase)
            verifyError(testCase,@()tsphopfieldnetOptions('Names',1:10),'MATLAB:InputParser:ArgumentFailedValidation');           
        end
        function tsphopfieldnetOptions_NamesNotCellStr_Errors(testCase)
            verifyError(testCase,@()tsphopfieldnetOptions('Names',num2cell(1:10)),'MATLAB:InputParser:ArgumentFailedValidation');           
        end
        function tsphopfieldnetOptions_DistanceTypeNotChar_Errors(testCase)
            verifyError(testCase,@()tsphopfieldnetOptions('DistanceType',pi),'MATLAB:InputParser:ArgumentFailedValidation');
        end
        function tsphopfieldnetOptions_DistanceTypeNotMatchValidOptions_Errors(testCase)
            verifyError(testCase,@()tsphopfieldnetOptions('DistanceType','euclidean'),'MATLAB:HopfieldNetworkTSPOptions:unrecognizedStringChoice');
        end
        function tsphopfieldnetOptions_DistTypeNotExplDistGivenCoordEmpty_Errors(testCase)
            verifyError(testCase,@()tsphopfieldnetOptions('DistanceType','euc','DistanceMatrix',rand(4),'Coordinates',[]),'HopfieldNetworkTSPOptions:DistanceMatrixGivenTypeNotExplicit');
        end
        function tsphopfieldnetOptions_DistTypeExplDistGivenCoordGiven_Errors(testCase)
            verifyError(testCase,@()tsphopfieldnetOptions('DistanceType','explicit','DistanceMatrix',rand(4),'Coordinates',rand(4,2)),'HopfieldNetworkTSPOptions:nonEmptyCoordinates');
        end
        function tsphopfieldnetOptions_DistTypeExplDistEmptyCoordGiven_Errors(testCase)
            verifyError(testCase,@()tsphopfieldnetOptions('DistanceType','explicit','DistanceMatrix',[],'Coordinates',rand(4,2)),'HopfieldNetworkTSPOptions:emptyDistanceMatrix');
        end
        function tsphopfieldnetOptions_CoordinatesNotReal_Errors(testCase)
            verifyError(testCase,@()tsphopfieldnetOptions('Coordinates',[1,2;-1,1;0,1+1i]),'MATLAB:InputParser:ArgumentFailedValidation');
        end
        function tsphopfieldnetOptions_CoordinatesNotNumeric_Errors(testCase)
            verifyError(testCase,@()tsphopfieldnetOptions('Coordinates','coords'),'MATLAB:InputParser:ArgumentFailedValidation');
        end
        function tsphopfieldnetOptions_CoordinatesNotTwoColumns_Errors(testCase)
            verifyError(testCase,@()tsphopfieldnetOptions('Coordinates',rand(2,10)),'MATLAB:InputParser:ArgumentFailedValidation');
        end
        function tsphopfieldnetOptions_SubtoursNotCell_Errors(testCase)
            verifyError(testCase,@()tsphopfieldnetOptions('Subtours',1:10),'MATLAB:InputParser:ArgumentFailedValidation');           
        end
        function tsphopfieldnetOptions_SubtoursNotCellStr_Errors(testCase)
            verifyError(testCase,@()tsphopfieldnetOptions('Subtours',num2cell(1:10)),'MATLAB:InputParser:ArgumentFailedValidation');           
        end
        function tsphopfieldnetOptions_SubtoursNotMatchSubtoursPositions_Errors(testCase)
            verifyError(testCase,@()tsphopfieldnetOptions('Subtours',{'A','B','C'}),'HopfieldNetworkTSPOptions:unevensubtours');
        end        
        function tsphopfieldnetOptions_SubtoursPositionsNotReal_Errors(testCase)
            verifyError(testCase,@()tsphopfieldnetOptions('SubtoursPositions',1+2i),'MATLAB:InputParser:ArgumentFailedValidation');
        end
        function tsphopfieldnetOptions_SubtoursPositionsNotNumeric_Errors(testCase)
            verifyError(testCase,@()tsphopfieldnetOptions('SubtoursPositions','1'),'MATLAB:InputParser:ArgumentFailedValidation');
        end
        function tsphopfieldnetOptions_SubtoursPositionsNotInteger_Errors(testCase)
            verifyError(testCase,@()tsphopfieldnetOptions('SubtoursPositions',[1,2,3.01]),'MATLAB:InputParser:ArgumentFailedValidation');           
        end
        function tsphopfieldnetOptions_SubtoursPositionsNotMatchSubtours_Errors(testCase)
            verifyError(testCase,@()tsphopfieldnetOptions('SubtoursPositions',[1,2,3]),'HopfieldNetworkTSPOptions:unevensubtours');
        end
        function tsphopfieldnetOptions_TauNotScalar_Errors(testCase)
            verifyError(testCase,@()tsphopfieldnetOptions('Tau',[10,10]),'MATLAB:InputParser:ArgumentFailedValidation');
        end
        function tsphopfieldnetOptions_TauNotReal_Errors(testCase)
            verifyError(testCase,@()tsphopfieldnetOptions('Tau',1+2i),'MATLAB:InputParser:ArgumentFailedValidation');           
        end
        function tsphopfieldnetOptions_TauNotNumeric_Errors(testCase)
            verifyError(testCase,@()tsphopfieldnetOptions('Tau','5'),'MATLAB:InputParser:ArgumentFailedValidation');           
        end       
        function tsphopfieldnetOptions_TauNotInteger_Errors(testCase)
            verifyError(testCase,@()tsphopfieldnetOptions('Tau',10.1),'MATLAB:InputParser:ArgumentFailedValidation');           
        end
        function tsphopfieldnetOptions_TauNotGreaterOrEqualThanZero_Errors(testCase)
            verifyError(testCase,@()tsphopfieldnetOptions('Tau',-1),'MATLAB:InputParser:ArgumentFailedValidation');           
        end
        function tsphopfieldnetOptions_UnmatchedParameter_Errors(testCase)
            verifyError(testCase,@()tsphopfieldnetOptions('InvTransferFcn', 'atanh'),'MATLAB:InputParser:UnmatchedParameter');
        end        
    end
end
