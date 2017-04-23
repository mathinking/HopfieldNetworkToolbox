classdef test_tsphopfieldnet < matlab.unittest.TestCase
%TEST_TSPHOPFIELDNET Summary of this class goes here
%   Detailed explanation goes here
    
	methods (Test)
        % [the name of the tested method]_[expected input / tested state]_[expected behavior]
        
        % Verify inputs
        % Checking tsphopfieldnet input arguments
        function tsphopfieldnet_InstanceDerivesHopfieldnetwork_Ok(testCase)
            verifyInstanceOf(testCase, tsphopfieldnet(4,1), ?network.HopfieldNetwork);
        end
        function tsphopfieldnet_InstanceIsTsphopfieldnet_Ok(testCase)
            verifyInstanceOf(testCase, tsphopfieldnet(4,1), ?network.HopfieldNetworkTSP);
        end   
        function tsphopfieldnet_IncorrecNumberInputs_Errors(testCase)
            verifyError(testCase,@()tsphopfieldnet(),'HopfieldNetworkTSP:IncorrectInputArguments');       
        end
        function tsphopfieldnet_networkSizeNotReal_Errors(testCase)
            verifyError(testCase,@()tsphopfieldnet(1+2i,0.1),'HopfieldNetwork:invalid_value');
        end
        function tsphopfieldnet_networkSizeNotNumeric_Errors(testCase)
            verifyError(testCase,@()tsphopfieldnet('size',0.1),'HopfieldNetwork:invalid_value');
        end
        function tsphopfieldnet_networkSizeNotScalar_Errors(testCase)
            verifyError(testCase,@()tsphopfieldnet(1:5,0.1),'HopfieldNetwork:invalid_value');
        end
        function tsphopfieldnet_networkSizeNotInteger_Errors(testCase)
            verifyError(testCase,@()tsphopfieldnet(10.2,0.1),'HopfieldNetwork:invalid_value');
        end
        function tsphopfieldnet_networkSizeNotGreaterOrEqualThanOne_Errors(testCase)
            verifyError(testCase,@()tsphopfieldnet(0,0.1),'HopfieldNetwork:invalid_value');
        end
        function tsphopfieldnet_networkParameterNotFloat_Errors(testCase)
            verifyError(testCase,@()tsphopfieldnet(10,int8(1)),'HopfieldNetworkTSP:invalid_value');
        end
        function tsphopfieldnet_networkParameterNotReal_Errors(testCase)
            verifyError(testCase,@()tsphopfieldnet(10,1+2i),'HopfieldNetworkTSP:invalid_value');
        end
        function tsphopfieldnet_networkParameterNotNumeric_Errors(testCase)
            verifyError(testCase,@()tsphopfieldnet(10,'C'),'HopfieldNetworkTSP:invalid_value');
        end
        function tsphopfieldnet_networkParameterNotScalar_Errors(testCase)
            verifyError(testCase,@()tsphopfieldnet(10,0:0.1:1),'HopfieldNetworkTSP:invalid_value');
        end
        function tsphopfieldnet_networkParameterNotGreaterThanZero_Errors(testCase)
            verifyError(testCase,@()tsphopfieldnet(10,0),'HopfieldNetworkTSP:invalid_value');
        end
        function tsphopfieldnet_optionsNotTsphopfieldnetOptions_Errors(testCase)
            verifyError(testCase,@()tsphopfieldnet(10,0.1,'options'),'HopfieldNetworkTSP:invalid_datatype');
        end

        % Verify options set
        function tsphopfieldnet_KNotLessThanFloorNDividedByTwo_Errors(testCase)
            N = 5;
            options = tsphopfieldnetOptions('K', 3);
            verifyError(testCase,@()tsphopfieldnet(N,0.1,options),'HopfieldNetworkTSP:invalid_value');
        end
        function tsphopfieldnet_NamesLengthNotMatchNetworkSize_Errors(testCase)
            N = 4;
            options = tsphopfieldnetOptions('Names', {'Madrid','London','New York'});
            verifyError(testCase,@()tsphopfieldnet(N,0.1,options),'HopfieldNetworkTSP:invalid_value');
        end
        function tsphopfieldnet_CoordinatesLengthNotMatchNetworkSize_Errors(testCase)
            coords = rand(5,2);
            options = tsphopfieldnetOptions('Coordinates', coords);
            verifyError(testCase,@()tsphopfieldnet(4,0.1,options),'HopfieldNetworkTSP:invalid_value');
        end
        function tsphopfieldnet_NoCoordsDistMatrixSizeNotMatchNetworkSize_Errors(testCase)
            d = rand(5);
            options = tsphopfieldnetOptions('DistanceMatrix', d, 'DistanceType', 'explicit');
            verifyError(testCase,@()tsphopfieldnet(4,0.1,options),'HopfieldNetworkTSP:invalid_value');
        end
        function tsphopfieldnet_CoordsDistMatrixSizeNotMatchNetworkSize_Errors(testCase)
            d = rand(5);
            options = tsphopfieldnetOptions('DistanceMatrix', d, 'DistanceType', 'euc', 'Coordinates', rand(5,2));
            verifyError(testCase,@()tsphopfieldnet(4,0.1,options),'HopfieldNetworkTSP:invalid_value');
        end
        function tsphopfieldnet_CoordinatesNonFloatDistMatrixEmpty_Warning(testCase)
            options = tsphopfieldnetOptions('Coordinates',int8([1,2;3,4;5,6;7,8;9,10;-1,-2;-3,-4;-5,-6;-7,-8;-9,-10]));
            verifyWarning(testCase, @()tsphopfieldnet(10,0.1,options),'HopfieldNetworkTSP:setOptions:DataConversion');
        end
        function tsphopfieldnet_CoordinatesNonFloatDistMatrixNotEmpty_Warning(testCase)
            options = tsphopfieldnetOptions('Coordinates',int8([1,2;3,4;5,6;7,8;9,10;-1,-2;-3,-4;-5,-6;-7,-8;-9,-10]),'DistanceMatrix',rand(10));
            verifyWarning(testCase, @()tsphopfieldnet(10,0.1,options),'HopfieldNetworkTSP:setOptions:DataConversion');
        end
        function tsphopfieldnet_DistanceMatrixNonFloatCoordsEmpty_Warning(testCase)
            options = tsphopfieldnetOptions('Coordinates',[],'DistanceMatrix',randi(100,10,10,'uint8'),'DistanceType','explicit');
            verifyWarning(testCase, @()tsphopfieldnet(10,0.1,options),'HopfieldNetworkTSP:setOptions:DataConversion');
        end
        function tsphopfieldnet_DistanceMatrixNonFloatCoordsNonEmpty_Warning(testCase)
            options = tsphopfieldnetOptions('Coordinates',rand(10,2),'DistanceMatrix',randi(100,10,10,'uint8'));
            verifyWarning(testCase, @()tsphopfieldnet(10,0.1,options),'HopfieldNetworkTSP:setOptions:DataConversion');
        end      
        function tsphopfieldnet_SubtourPositionGreaterThanNetworkSize_Errors(testCase)
            options = tsphopfieldnetOptions('Subtours', {'A-B-C'}, 'SubtoursPositions', 5);
            verifyError(testCase,@()tsphopfieldnet(4,0.1,options),'HopfieldNetworkTSP:invalidsubtoursPositions');
        end
        function tsphopfieldnet_SubtourCityNotPresentInNames_Errors(testCase)
            options = tsphopfieldnetOptions('Names',{'A','B','C','D','E'},'Subtours', {'H-C'}, 'SubtoursPositions', 5);
            verifyError(testCase,@()tsphopfieldnet(5,0.1,options),'HopfieldNetworkTSP:SubtourCityNotPresentInNames');
        end
        function tsphopfieldnet_SubtourCityRepeatedSameSubtour_Errors(testCase)
            options = tsphopfieldnetOptions('Subtours', {'A-B-A'}, 'SubtoursPositions', 1);
            verifyError(testCase,@()tsphopfieldnet(5,0.1,options),'HopfieldNetworkTSP:SubtourCityRepeatedSameSubtour');
        end
        function tsphopfieldnet_SubtourCitiesRepeatedSameSubtour_Errors(testCase)
            options = tsphopfieldnetOptions('Subtours', {'A-B-A-B'}, 'SubtoursPositions', 1);
            verifyError(testCase,@()tsphopfieldnet(5,0.1,options),'HopfieldNetworkTSP:SubtourCitiesRepeatedSameSubtour');
        end
        function tsphopfieldnet_SubtourCityRepeatedInDifferentSubtours_Errors(testCase)
            options = tsphopfieldnetOptions('Subtours', {'A-B-C', 'B-D'}, 'SubtoursPositions', [1,4]);
            verifyError(testCase,@()tsphopfieldnet(4,0.1,options),'HopfieldNetworkTSP:verifyIfValidSubtours:repeatedCity');
        end
        function tsphopfieldnet_RepeatedPositionInSubtours_Errors(testCase)
            options = tsphopfieldnetOptions('Subtours', {'A-B-C', 'D-E'}, 'SubtoursPositions', [1,3]);
            verifyError(testCase,@()tsphopfieldnet(5,0.1,options),'HopfieldNetworkTSP:verifyIfValidSubtours:repeatedSubtourPosition');
        end
        function tsphopfieldnet_TauGreaterThanNminusOne_Errors(testCase)
            N = 10;
            options = tsphopfieldnetOptions('Tau',N);
            verifyError(testCase,@()tsphopfieldnet(N,0.1,options),'HopfieldNetworkTSP:invalid_value');
        end

        % Checking f^-1(f(x)) = x, f=satlin
        function tsphopfieldnet_TransferFcnIssatlin_InvTransferFcnIsInvsatlin(testCase)
            options = tsphopfieldnetOptions('U0',0.5,'TransferFcn','satlin');
            networkSize = 10;
            C = 0.1;           
            net = tsphopfieldnet(networkSize, C, options);
            transferFcn = getSetting(net,'TransferFcn');
            invTransferFcn = getSetting(net,'InvTransferFcn');
            value = 0.1234;
            verifyEqual(testCase, invTransferFcn(transferFcn(value)), value, 'AbsTol', power(10, -1 * getSetting(net,'E')));
        end
        % Checking f^-1(f(x)) = x, f=tanh        
        function tsphopfieldnet_transferFcnIstanh_invTransferFcnIsatanh(testCase)
            options = tsphopfieldnetOptions('U0',0.5,'TransferFcn','tanh');
            networkSize = 10;
            C = 0.1;           
            net = tsphopfieldnet(networkSize, C, options);
            transferFcn = getSetting(net,'TransferFcn');
            invTransferFcn = getSetting(net,'InvTransferFcn');
            value = 0.1234;
            verifyEqual(testCase, invTransferFcn(transferFcn(value)), value, 'AbsTol', power(10, -1 * getSetting(net,'E')));
        end

%         function createOptions_coordsMatrixNotDouble_Errors(testCase)
%             verifyError(testCase, @()tsphopfieldnet.createOptions('coords',rand(4,2,'single')),'HopfieldNetworkTSP:invalid_datatype');       
%         end
%         function createOptions_coordsMatrixNotTwoColumns_Errors(testCase)
%             verifyError(testCase, @()tsphopfieldnet.createOptions('coords',rand(4,3)),'HopfieldNetworkTSP:invalid_value');       
%         end
%         function setOptions_coordsMatrixNotNetworkSize_Errors(testCase)
%             options = tsphopfieldnet.createOptions('coords',rand(4,2));
%             networkSize = 10;
%             C = 0.1;
%             verifyError(testCase, @()tsphopfieldnet(networkSize, C, options),'HopfieldNetworkTSP:invalid_value');       
%         end
%         function setOptions_cityNamesNotCell_Errors(testCase)
%             verifyError(testCase, @()tsphopfieldnet.createOptions('names','Madrid'),'HopfieldNetworkTSP:invalid_datatype');       
%         end
%         function setOptions_cityNamesNotCellArrayOfChars_Errors(testCase)
%             options = tsphopfieldnet.createOptions('names',{1,'Madrid','2'});
%             networkSize = 3;
%             C = 0.1;
%             verifyError(testCase, @()tsphopfieldnet(networkSize, C, options),'HopfieldNetworkTSP:invalid_value');       
%         end
%         function setOptions_cityNamesNotNetworkSize_Errors(testCase)
%             options = tsphopfieldnet.createOptions('names',{'London','Madrid','Berlin','Paris'});
%             networkSize = 3;
%             C = 0.1;
%             verifyError(testCase, @()tsphopfieldnet(networkSize, C, options),'HopfieldNetworkTSP:invalid_value');       
%         end   
%         
%         
%         function createOptions_trainFcnNotChar_Errors(testCase)
%             verifyError(testCase, @()tsphopfieldnet.createOptions('trainFcn',{'trainty'}),'HopfieldNetworkTSP:invalid_datatype');          
%         end
%         function createOptions_trainFcnNotValid_Errors(testCase)
%             verifyError(testCase, @()tsphopfieldnet.createOptions('trainFcn','trainyt'),'HopfieldNetworkTSP:invalid_value');          
%         end
%         function createOptions_simFcnNotChar_Errors(testCase)
%             verifyError(testCase, @()tsphopfieldnet.createOptions('simFcn',{'talavan-yanez'}),'HopfieldNetworkTSP:invalid_datatype');          
%         end
%         function createOptions_simFcnNotValid_Errors(testCase)
%             verifyError(testCase, @()tsphopfieldnet.createOptions('simFcn','yanez-talavan'),'HopfieldNetworkTSP:invalid_value');          
%         end
%         function createOptions_simFcnTalavanYanez(testCase)
%         	options = tsphopfieldnet.createOptions('simFcn','talavan-yanez');
%             net = tsphopfieldnet(4, 1e-5, options);
%             verifyTrue(testCase, strcmp(getSimFcn(net),options.simFcn));
%         end
%         function createOptions_simFcnDivideConquer(testCase)
%         	options = tsphopfieldnet.createOptions('simFcn','divide-conquer');
%             net = tsphopfieldnet(4, 1e-5, options);
%             verifyTrue(testCase, strcmp(getSimFcn(net),options.simFcn));
%         end        
%         
%         function train_trainFcnIsTrainty_WorksFine(testCase)
%             import matlab.unittest.constraints.IsEqualTo;
%             import matlab.unittest.constraints.AbsoluteTolerance;
%                         
%             networkSize = 6;
%             C = 0.1;
%             options = tsphopfieldnet.createOptions('trainFcn','trainty');
%             net = tsphopfieldnet(networkSize, C, options);
%             trainParamExpected = getTrainParam(net);
%             trainParamExpected.A = 2.6;
%             trainParamExpected.B = 3.1;
%             trainParamExpected.D = 1;
%             trainParamExpected.Np = 36;
%             trainParamExpected.dL = 0.5;
%             trainParamExpected.dU = 1;
%             trainParamExpected.dUaux = 2;
%             trainParamExpected.K = 0;
%             trainParamExpected.rho = 0.5;
%             trainParamExpected = orderfields(trainParamExpected);
%             train(net);
%             trainParam = getTrainParam(net);
%             verifyThat(testCase, trainParam, IsEqualTo(trainParamExpected, 'Within', AbsoluteTolerance(power(10, -1 * getSetting(net,'e')))));
%         end
%         function train_trainFcnIsNotTrainty_Errors(testCase)
%             networkSize = 6;
%             C = 0.1;
%             options = tsphopfieldnet.createOptions('trainFcn','trainty');
%             options.trainFcn = 'nottrainty';
%             net = tsphopfieldnet(networkSize, C, options);
%             verifyError(testCase, @()train(net),'HopfieldNetworkTSP:unvalidTrainFcn');
%         end
%         
%         function sim_simTalavanYanezOutputWithPolygon_hasOne1perRow(testCase)
%             rng(2);
%             options = tsphopfieldnet.createOptions('simFcn','talavan-yanez');
%             networkSize = 6;
%             C = 0.1;
%             net = tsphopfieldnet(networkSize, C, options);
%             train(net);
%             V = sim(net);
%             verifyTrue(testCase, all(sum(V,2) == 1));
%         end
%         function sim_simTalavanYanezOutputWithPolygon_hasOne1perCol(testCase)
%             rng(2);
%             options = tsphopfieldnet.createOptions('simFcn','talavan-yanez');
%             networkSize = 6;
%             C = 0.1;
%             net = tsphopfieldnet(networkSize, C, options);
%             train(net);
%             V = sim(net);
%             verifyTrue(testCase, all(sum(V,1) == 1));
%         end
%         function sim_simTalavanYanezOutputWithPolygon_sumsN(testCase)
%             rng(2);
%             options = tsphopfieldnet.createOptions('simFcn','talavan-yanez');
%             networkSize = 6;
%             C = 0.1;
%             net = tsphopfieldnet(networkSize, C, options);
%             train(net);
%             V = sim(net);
%             verifyEqual(testCase, sum(sum(V)), networkSize);
%         end        
% 
%         function sim_simDivideConquerOutputWithPolygon_hasOne1perRow(testCase)
%             rng(2);
%             options = tsphopfieldnet.createOptions('simFcn','divide-conquer');
%             networkSize = 6;
%             C = 0.1;
%             net = tsphopfieldnet(networkSize, C, options);
%             train(net);
%             V = sim(net);
%             verifyTrue(testCase, all(sum(V,2) == 1));
%         end        
%         function sim_simDivideConquerOutputWithPolygon_hasOne1perCol(testCase)
%             rng(2);
%             options = tsphopfieldnet.createOptions('simFcn','divide-conquer');
%             networkSize = 6;
%             C = 0.1;
%             net = tsphopfieldnet(networkSize, C, options);
%             train(net);
%             V = sim(net);
%             verifyTrue(testCase, all(sum(V,1) == 1));
%         end
%         function sim_simDivideConquerOutputWithPolygon_sumsN(testCase)
%             rng(2);
%             options = tsphopfieldnet.createOptions('simFcn','divide-conquer');
%             networkSize = 6;
%             C = 0.1;
%             net = tsphopfieldnet(networkSize, C, options);
%             train(net);
%             V = sim(net);
%             verifyEqual(testCase, sum(sum(V)), networkSize);
%         end    
%         
%         function sim_simTalavanYanezOutputWithBerlin52_hasOne1perRow(testCase)
%             rng(2);
%             problem = tsplib({'berlin52'});
%             options = tsphopfieldnet.createOptions('coords',problem.coords,'type',problem.type,'simFcn','talavan-yanez');
%             networkSize = problem.nCities;
%             C = 0.1;
%             net = tsphopfieldnet(networkSize, C, options);
%             train(net);
%             V = sim(net);
%             verifyTrue(testCase, all(sum(V,2) == 1));
%         end
%         function sim_simTalavanYanezOutputWithBerlin52_hasOne1perCol(testCase)
%             rng(2);
%             problem = tsplib({'berlin52'});
%             options = tsphopfieldnet.createOptions('coords',problem.coords,'type',problem.type,'simFcn','talavan-yanez');
%             networkSize = problem.nCities;
%             C = 0.1;
%             net = tsphopfieldnet(networkSize, C, options);
%             train(net);
%             V = sim(net);
%             verifyTrue(testCase, all(sum(V,1) == 1));
%         end
%         function sim_simTalavanYanezOutputWithBerlin52_sumsN(testCase)
%             rng(2);
%             problem = tsplib({'berlin52'});
%             options = tsphopfieldnet.createOptions('coords',problem.coords,'type',problem.type,'simFcn','talavan-yanez');
%             networkSize = problem.nCities;
%             C = 0.1;
%             net = tsphopfieldnet(networkSize, C, options);
%             train(net);
%             V = sim(net);
%             verifyEqual(testCase, sum(sum(V)), networkSize);
%         end
%         
%         function sim_simDivideConquerOutputWithBerlin52_hasOne1perRow(testCase)
%             rng(2);
%             problem = tsplib({'berlin52'});
%             options = tsphopfieldnet.createOptions('coords',problem.coords,'type',problem.type,'simFcn','divide-conquer');
%             networkSize = problem.nCities;
%             C = 0.1;
%             net = tsphopfieldnet(networkSize, C, options);
%             train(net);
%             V = sim(net);
%             verifyTrue(testCase, all(sum(V,2) == 1));
%         end
%         function sim_simDivideConquerOutputWithBerlin52_hasOne1perCol(testCase)
%             rng(2);
%             problem = tsplib({'berlin52'});
%             options = tsphopfieldnet.createOptions('coords',problem.coords,'type',problem.type,'simFcn','divide-conquer');
%             networkSize = problem.nCities;
%             C = 0.1;
%             net = tsphopfieldnet(networkSize, C, options);
%             train(net);
%             V = sim(net);
%             verifyTrue(testCase, all(sum(V,1) == 1));
%         end
%         function sim_simDivideConquerOutputWithBerlin52_sumsN(testCase)
%             rng(2);
%             problem = tsplib({'berlin52'});
%             options = tsphopfieldnet.createOptions('coords',problem.coords,'type',problem.type,'simFcn','divide-conquer');
%             networkSize = problem.nCities;
%             C = 0.1;
%             net = tsphopfieldnet(networkSize, C, options);
%             train(net);
%             V = sim(net);
%             verifyEqual(testCase, sum(sum(V)), networkSize);
%         end             
%         
%         function reinit_netObjectAlreadySimulated_WorksFine(testCase)
%             import matlab.unittest.constraints.IsEqualTo;
%                         
%             networkSize = 6;
%             C = 0.1;
%             net = tsphopfieldnet(networkSize, C);
%             expectedResults = getResults(net);
%             train(net);
%             sim(net);
%             reinit(net);
%             results = getResults(net);
%             verifyThat(testCase, results, IsEqualTo(expectedResults));
%         end
%         
%         function saddle_simFcnNotTalavanYanez_saddleGetsComputed(testCase)
%             networkSize = 10;
%             C = 0.0001;
%             options = tsphopfieldnet.createOptions('simFcn', 'divide-conquer');
%             net = tsphopfieldnet(networkSize, C, options);
%             simFcnExpected = getSimFcn(net);
%             saddle(net);
%             verifyEqual(testCase, getSimFcn(net), simFcnExpected);
%         end
%         function saddle_computingSaddle_allColumnsEqual(testCase)
%             import matlab.unittest.constraints.IsEqualTo;
%             import matlab.unittest.constraints.AbsoluteTolerance;
% 
%             problem = tsplib({'berlin52'});
%             options = tsphopfieldnet.createOptions('coords',problem.coords,'type',problem.type,'simFcn','talavan-yanez');
%             networkSize = problem.nCities;
%             C = 0.1;
%             net = tsphopfieldnet(networkSize, C, options);
%             train(net);
%             V = saddle(net);
%             for i = 2:networkSize
%                 verifyThat(testCase, V(:,1), IsEqualTo(V(:,i), 'Within', AbsoluteTolerance(power(10, -1 * getSetting(net,'e')))));
%             end
%         end
        
     end
 end
