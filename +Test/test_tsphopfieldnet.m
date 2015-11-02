classdef test_tsphopfieldnet < matlab.unittest.TestCase
%TEST_TSPHOPFIELDNET Summary of this class goes here
%   Detailed explanation goes here
    
	methods (Test)
        % [the name of the tested method]_[expected input / tested state]_[expected behavior]
        
        % Verify inputs
        % Checking tsphopfieldnet input arguments
        function tsphopfieldnet_IncorrecNumberInputs_Errors(testCase)
            verifyError(testCase,@()tsphopfieldnet(),'tsphopfieldnet:IncorrectInputArguments');       
        end
        function tsphopfieldnet_networkSizeNotDouble_Errors(testCase)
            verifyError(testCase,@()tsphopfieldnet(uint8(10),0.1),'hopfieldnetwork:invalid_datatype');
        end
        function tsphopfieldnet_networkSizeTooSmall_Errors(testCase)
            verifyError(testCase,@()tsphopfieldnet(1,0.1),'tsphopfieldnet:invalid_value');
        end
        function tsphopfieldnet_CNotDouble_Errors(testCase)
            verifyError(testCase,@()tsphopfieldnet(10,single(0.1)),'tsphopfieldnet:invalid_datatype');
        end
        function tsphopfieldnet_CNotPositive_Errors(testCase)
            verifyError(testCase,@()tsphopfieldnet(10,0),'tsphopfieldnet:invalid_value');
        end
        function tsphopfieldnet_optionsNotStruct_Errors(testCase)
            verifyError(testCase,@()tsphopfieldnet(10,0.1,'options'),'hopfieldnetwork:invalid_datatype');
        end        
       
        function createOptions_transferFcnIssatlin_invTransferFcnIsinvsatlin(testCase)
            options = tsphopfieldnet.createOptions('u0',0.5,'transferFcn','satlin');
            networkSize = 10;
            C = 0.1;           
            net = tsphopfieldnet(networkSize, C, options);
            transferFcn = getSetting(net,'transferFcn');
            invTransferFcn = getSetting(net,'invTransferFcn');
            value = 0.1234;
            verifyEqual(testCase, invTransferFcn(transferFcn(value)), value, 'AbsTol', power(10, -1 * getSetting(net,'e')));
        end
        
        function createOptions_invTransferFcnIsinvsatlin_transferFcnIssatlin(testCase)
            options = tsphopfieldnet.createOptions('u0',0.5,'invTransferFcn','invsatlin');
            networkSize = 10;
            C = 0.1;           
            net = tsphopfieldnet(networkSize, C, options);
            invTransferFcn = getSetting(net,'invTransferFcn');
            transferFcn = getSetting(net,'transferFcn');
            value = 0.1234;
            verifyEqual(testCase, transferFcn(invTransferFcn(value)), value, 'AbsTol', power(10, -1 * getSetting(net,'e')));
        end
        
        function createOptions_transferFcnIstanh_invTransferFcnIsatanh(testCase)
            options = tsphopfieldnet.createOptions('u0',0.5,'transferFcn','tanh');
            networkSize = 10;
            C = 0.1;           
            net = tsphopfieldnet(networkSize, C, options);
            transferFcn = getSetting(net,'transferFcn');
            invTransferFcn = getSetting(net,'invTransferFcn');
            value = 0.1234;
            verifyEqual(testCase, invTransferFcn(transferFcn(value)), value, 'AbsTol', power(10, -1 * getSetting(net,'e')));
        end

        function createOptions_invTransferFcnIsatanh_transferFcnIstanh(testCase)
            options = tsphopfieldnet.createOptions('u0',0.5,'invTransferFcn','atanh');
            networkSize = 10;
            C = 0.1;           
            net = tsphopfieldnet(networkSize, C, options);
            invTransferFcn = getSetting(net,'invTransferFcn');
            transferFcn = getSetting(net,'transferFcn');
            value = 0.1234;
            verifyEqual(testCase, transferFcn(invTransferFcn(value)), value, 'AbsTol', power(10, -1 * getSetting(net,'e')));
        end
         
        function createOptions_distanceMatrixNotDouble_Errors(testCase)
            verifyError(testCase, @()tsphopfieldnet.createOptions('d',rand(4,'single')),'tsphopfieldnet:invalid_datatype');       
        end        
        function createOptions_distanceMatrixNotSquare_Errors(testCase)
            verifyError(testCase, @()tsphopfieldnet.createOptions('d',rand(3,4)),'tsphopfieldnet:invalid_value');       
        end
        function tsphopfieldnet_distanceMatrixNotSameSizeNetwork_Errors(testCase)
            options = tsphopfieldnet.createOptions('d',rand(4));
            networkSize = 10;
            C = 0.1;
            verifyError(testCase, @()tsphopfieldnet(networkSize, C, options),'tsphopfieldnet:invalid_value');       
        end

        function createOptions_coordsMatrixNotDouble_Errors(testCase)
            verifyError(testCase, @()tsphopfieldnet.createOptions('coords',rand(4,2,'single')),'tsphopfieldnet:invalid_datatype');       
        end
        function createOptions_coordsMatrixNotTwoColumns_Errors(testCase)
            verifyError(testCase, @()tsphopfieldnet.createOptions('coords',rand(4,3)),'tsphopfieldnet:invalid_value');       
        end
        function setOptions_coordsMatrixNotNetworkSize_Errors(testCase)
            options = tsphopfieldnet.createOptions('coords',rand(4,2));
            networkSize = 10;
            C = 0.1;
            verifyError(testCase, @()tsphopfieldnet(networkSize, C, options),'tsphopfieldnet:invalid_value');       
        end
        function setOptions_cityNamesNotCell_Errors(testCase)
            verifyError(testCase, @()tsphopfieldnet.createOptions('names','Madrid'),'tsphopfieldnet:invalid_datatype');       
        end
        function setOptions_cityNamesNotCellArrayOfChars_Errors(testCase)
            options = tsphopfieldnet.createOptions('names',{1,'Madrid','2'});
            networkSize = 3;
            C = 0.1;
            verifyError(testCase, @()tsphopfieldnet(networkSize, C, options),'tsphopfieldnet:invalid_value');       
        end
        function setOptions_cityNamesNotNetworkSize_Errors(testCase)
            options = tsphopfieldnet.createOptions('names',{'London','Madrid','Berlin','Paris'});
            networkSize = 3;
            C = 0.1;
            verifyError(testCase, @()tsphopfieldnet(networkSize, C, options),'tsphopfieldnet:invalid_value');       
        end   
        
        % TODO tests to ensure that fixedCities and startinPos have the same length
        
        function createOptions_trainFcnNotChar_Errors(testCase)
            verifyError(testCase, @()tsphopfieldnet.createOptions('trainFcn',{'trainty'}),'tsphopfieldnet:invalid_datatype');          
        end
        function createOptions_trainFcnNotValid_Errors(testCase)
            verifyError(testCase, @()tsphopfieldnet.createOptions('trainFcn','trainyt'),'tsphopfieldnet:invalid_value');          
        end
        function createOptions_simFcnNotChar_Errors(testCase)
            verifyError(testCase, @()tsphopfieldnet.createOptions('simFcn',{'talavan-yanez'}),'tsphopfieldnet:invalid_datatype');          
        end
        function createOptions_simFcnNotValid_Errors(testCase)
            verifyError(testCase, @()tsphopfieldnet.createOptions('simFcn','yanez-talavan'),'tsphopfieldnet:invalid_value');          
        end
        
        function train_trainFcnIsTrainty_WorksFine(testCase)
            import matlab.unittest.constraints.IsEqualTo;
            import matlab.unittest.constraints.AbsoluteTolerance;
                        
            networkSize = 6;
            C = 0.1;
            options = tsphopfieldnet.createOptions('trainFcn','trainty');
            net = tsphopfieldnet(networkSize, C, options);
            trainParamExpected = getTrainParam(net);
            trainParamExpected.A = 2.6;
            trainParamExpected.B = 3.1;
            trainParamExpected.D = 1;
            trainParamExpected.Np = 36;
            trainParamExpected.dL = 0.5;
            trainParamExpected.dU = 1;
            trainParamExpected.dUaux = 2;
            trainParamExpected.K = 0;
            trainParamExpected.rho = 0.5;
            trainParamExpected = orderfields(trainParamExpected);
            train(net);
            trainParam = getTrainParam(net);
            verifyThat(testCase, trainParam, IsEqualTo(trainParamExpected, 'Within', AbsoluteTolerance(power(10, -1 * getSetting(net,'e')))));
        end
        function train_trainFcnIsNotTrainty_Errors(testCase)
            networkSize = 6;
            C = 0.1;
            options = tsphopfieldnet.createOptions('trainFcn','trainty');
            options.trainFcn = 'nottrainty';
            net = tsphopfieldnet(networkSize, C, options);
            verifyError(testCase, @()train(net),'tsphopfieldnet:unvalidTrainFcn');
        end
        
        function sim_simTalavanYanezOutputWithPolygon_hasOne1perRow(testCase)
            rng(2);
            options = tsphopfieldnet.createOptions('simFcn','talavan-yanez');
            networkSize = 6;
            C = 0.1;
            net = tsphopfieldnet(networkSize, C, options);
            train(net);
            V = sim(net);
            verifyTrue(testCase, all(sum(V,2) == 1));
        end
        function sim_simTalavanYanezOutputWithPolygon_hasOne1perCol(testCase)
            rng(2);
            options = tsphopfieldnet.createOptions('simFcn','talavan-yanez');
            networkSize = 6;
            C = 0.1;
            net = tsphopfieldnet(networkSize, C, options);
            train(net);
            V = sim(net);
            verifyTrue(testCase, all(sum(V,1) == 1));
        end
        function sim_simTalavanYanezOutputWithPolygon_sumsN(testCase)
            rng(2);
            options = tsphopfieldnet.createOptions('simFcn','talavan-yanez');
            networkSize = 6;
            C = 0.1;
            net = tsphopfieldnet(networkSize, C, options);
            train(net);
            V = sim(net);
            verifyEqual(testCase, sum(sum(V)), networkSize);
        end        

        function sim_simTalavanYanezOutputWithBerlin52_hasOne1perRow(testCase)
            rng(2);
            problem = tsplib({'berlin52'});
            options = tsphopfieldnet.createOptions('coords',problem.coords,'type',problem.type,'simFcn','talavan-yanez');
            networkSize = problem.nCities;
            C = 0.1;
            net = tsphopfieldnet(networkSize, C, options);
            train(net);
            V = sim(net);
            verifyTrue(testCase, all(sum(V,2) == 1));
        end
        function sim_simTalavanYanezOutputWithBerlin52_hasOne1perCol(testCase)
            rng(2);
            problem = tsplib({'berlin52'});
            options = tsphopfieldnet.createOptions('coords',problem.coords,'type',problem.type,'simFcn','talavan-yanez');
            networkSize = problem.nCities;
            C = 0.1;
            net = tsphopfieldnet(networkSize, C, options);
            train(net);
            V = sim(net);
            verifyTrue(testCase, all(sum(V,1) == 1));
        end
        function sim_simTalavanYanezOutputWithBerlin52_sumsN(testCase)
            rng(2);
            problem = tsplib({'berlin52'});
            options = tsphopfieldnet.createOptions('coords',problem.coords,'type',problem.type,'simFcn','talavan-yanez');
            networkSize = problem.nCities;
            C = 0.1;
            net = tsphopfieldnet(networkSize, C, options);
            train(net);
            V = sim(net);
            verifyEqual(testCase, sum(sum(V)), networkSize);
        end        
        
        function reinit_netObjectAlreadySimulated_WorksFine(testCase)
            import matlab.unittest.constraints.IsEqualTo;
                        
            networkSize = 6;
            C = 0.1;
            net = tsphopfieldnet(networkSize, C);
            expectedResults = getResults(net);
            train(net);
            sim(net);
            reinit(net);
            results = getResults(net);
            verifyThat(testCase, results, IsEqualTo(expectedResults));
        end
        
        function saddle_simFcnNotTalavanYanez_saddleGetsComputed(testCase)
            networkSize = 10;
            C = 0.0001;
            options = tsphopfieldnet.createOptions('simFcn', 'divide-conquer');
            net = tsphopfieldnet(networkSize, C, options);
            simFcnExpected = getSimFcn(net);
            saddle(net);
            verifyEqual(testCase, getSimFcn(net), simFcnExpected);
        end
        function saddle_computingSaddle_allColumnsEqual(testCase)
            import matlab.unittest.constraints.IsEqualTo;
            import matlab.unittest.constraints.AbsoluteTolerance;

            problem = tsplib({'berlin52'});
            options = tsphopfieldnet.createOptions('coords',problem.coords,'type',problem.type,'simFcn','talavan-yanez');
            networkSize = problem.nCities;
            C = 0.1;
            net = tsphopfieldnet(networkSize, C, options);
            train(net);
            V = saddle(net);
            for i = 2:networkSize
                verifyThat(testCase, V(:,1), IsEqualTo(V(:,i), 'Within', AbsoluteTolerance(power(10, -1 * getSetting(net,'e')))));
            end
        end
        
    end
end
