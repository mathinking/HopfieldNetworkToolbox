classdef test_hopfieldnetwork < matlab.unittest.TestCase
    %TESTTSPHOPFIELDNET Summary of this class goes here
    %   Detailed explanation goes here
    
	methods (Test)
        % [the name of the tested method]_[expected input / tested state]_[expected behavior]
                
        % Checking hopfieldnet method createOptions
%         function createOptions_Unkwown_Errors(testCase)
%              verifyError(testCase,@()hopfieldnetwork.createOptions('unknown',2),'hopfieldNetwork:unvalidSetting');       
%         end
%         
%         function createOptions_u0WrongType_Errors(testCase)
%              verifyError(testCase,@()hopfieldnetwork.createOptions('u0',uint8(2)),'hopfieldnetwork:invalid_datatype');
%         end
%         function createOptions_ExecutionEnvironmentWrongType_Errors(testCase)
%              verifyError(testCase,@()hopfieldnetwork.createOptions('ExecutionEnvironment',0),'hopfieldnetwork:invalid_datatype');
%         end
%         function createOptions_maxIterWrongType_Errors(testCase)
%              verifyError(testCase,@()hopfieldnetwork.createOptions('maxIter',true),'hopfieldnetwork:invalid_datatype');
%         end
%         function createOptions_epsWrongType_Errors(testCase)
%              verifyError(testCase,@()hopfieldnetwork.createOptions('e',{7}),'hopfieldnetwork:invalid_datatype');
%         end
%         function createOptions_qWrongType_Errors(testCase)
%              verifyError(testCase,@()hopfieldnetwork.createOptions('q',single(0.05)),'hopfieldnetwork:invalid_datatype');
%         end
%         function createOptions_R_ITERWrongType_Errors(testCase)
%              verifyError(testCase,@()hopfieldnetwork.createOptions('R_ITER',uint32(20)),'hopfieldnetwork:invalid_datatype');
%         end
%         function createOptions_dtWrongType_Errors(testCase)
%              verifyError(testCase,@()hopfieldnetwork.createOptions('dt',single(0.01)),'hopfieldnetwork:invalid_datatype');
%         end
%         function createOptions_verboseWrongType_Errors(testCase)
%              verifyError(testCase,@()hopfieldnetwork.createOptions('Verbose',1),'hopfieldnetwork:invalid_datatype');
%         end
%         function createOptions_transferFcnWrongType_Errors(testCase)
%              verifyError(testCase,@()hopfieldnetwork.createOptions('transferFcn',{'exponential'}),'hopfieldnetwork:invalid_datatype');
%         end
%         function createOptions_CheckpointPathWrongType_Errors(testCase)
%             verifyError(testCase,@()hopfieldnetwork.createOptions('CheckpointPath',1),'hopfieldnetwork:invalid_datatype');
%         end
%         function createOptions_SimulationPlotWrongType_Errors(testCase)
%             verifyError(testCase,@()hopfieldnetwork.createOptions('SimulationPlot',1),'hopfieldnetwork:invalid_datatype');
%         end
%         function createOptions_SimulationPlotPauseTimeWrongType_Errors(testCase)
%             verifyError(testCase,@()hopfieldnetwork.createOptions('SimulationPlotPauseTime',single(0.2)),'hopfieldnetwork:invalid_datatype');
%         end
% 
%         function createOptions_u0WrongValue_Errors(testCase)
%              verifyError(testCase,@()hopfieldnetwork.createOptions('u0',0),'hopfieldnetwork:invalid_value');
%         end
%         function createOptions_ExecutionEnvironmentWrongValue_Errors(testCase)
%              verifyError(testCase,@()hopfieldnetwork.createOptions('ExecutionEnvironment','FPGA'),'hopfieldnetwork:invalid_value');
%         end
%         function createOptions_maxIterWrongValue_Errors(testCase)
%              verifyError(testCase,@()hopfieldnetwork.createOptions('maxIter',0),'hopfieldnetwork:invalid_value');
%         end
%         function createOptions_eWrongValue_Errors(testCase)
%              verifyError(testCase,@()hopfieldnetwork.createOptions('e',50),'hopfieldnetwork:invalid_value');
%         end
%         function createOptions_qWrongValue_Errors(testCase)
%              verifyError(testCase,@()hopfieldnetwork.createOptions('q',0),'hopfieldnetwork:invalid_value');
%         end
%         function createOptions_R_ITERWrongValue_Errors(testCase)
%              verifyError(testCase,@()hopfieldnetwork.createOptions('R_ITER',-1),'hopfieldnetwork:invalid_value');
%         end
%         function createOptions_dtWrongValue_Errors(testCase)
%              verifyError(testCase,@()hopfieldnetwork.createOptions('dt',0),'hopfieldnetwork:invalid_value');
%         end        
%         function createOptions_transferFcnWrongValue_Errors(testCase)
%              verifyError(testCase,@()hopfieldnetwork.createOptions('transferFcn','linear'),'hopfieldnetwork:invalid_value');
%         end
%         function createOptions_SimulationPlotPauseTimeWrongValue_Errors(testCase)
%              verifyError(testCase,@()hopfieldnetwork.createOptions('SimulationPlotPauseTime',1),'hopfieldnetwork:invalid_value');
%         end
%         
%         function createOptions_ValidProperty_Works(testCase)
%             u0 = 0.35;
%             options = hopfieldnetwork.createOptions('u0',u0);
%             verifyEqual(testCase,options.setting.u0,u0);
%         end              
    end
    
end
