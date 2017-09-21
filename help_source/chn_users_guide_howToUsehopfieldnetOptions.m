%% How to use hopfieldnet and hopfieldnetOptions
% |hopfieldnet| is a function that returns an object of the class
% |HopfieldNetworkGQKP|. This object is a Continuous Hopfield Network
% designed to solve the Generalized Quadratic Knapsack Problem.
% 
% At the same time, |hopfieldnetOptions| returns an object of the class
% |HopfieldNetworkGQKPOptions|. Such object may be used to customize the
% Hopfield network object from |HopfieldNetworkGQKP|.
%
%% How |hopfieldnet| works
% |hopfieldnet|'s default behavior is to take 6 input arguments:
%%
% * P: Square matrix defining the quadratic term of the optimization problem
% * q: Vector defining the linear term of the optimization problem
% * Aeq: Matrix defining the linear equality constraints. Aeq*x = beq
% * beq: Vector defining the linear equality constraints. Aeq*x = beq
% * A: Matrix defining the linear inequality constraints. A*x <= b
% * b: Vector defining the linear inequality constraints. A*x <= b
% 
% Optionally, |hopfieldnet| may receive a seven input argument (of the 
% class |HopfieldNetworkGQKPOptions| using the function 
% |hopfieldnetOptions|) which gives further information of the GQKP
% problem considered. This way an arbitrary GQKP problem can be created. 
%
%%
% *Syntax*
%
%      net = hopfieldnet(P,q,Aeq,beq,A,b)
%      net = hopfieldnet(P,q,Aeq,beq,A,b,options)
% 
%% How |hopfieldnetOptions| works
% |hopfieldnetOptions| is a function returning an object of the class 
% |HopfieldNetworkGQKPOptions| which allows to specify all the settings and 
% details of the GQKP problem.
% 
%%
% *Syntax*
%
%      options = hopfieldnetOptions('param1','value1','param2','value2',...)

%%
% The following table lists the available options for 
% |HopfieldNetworkGQKPOptions|.
%
% <html>
% <style type="text/css">
% .tg  {border-collapse:collapse;border-spacing:0;}
% .tg td{font-family:Arial, sans-serif;font-size:14px;padding:10px 5px;border-style:solid;border-width:1px;overflow:hidden;word-break:normal;}
% .tg th{font-family:Arial, sans-serif;font-size:14px;font-weight:normal;padding:10px 5px;border-style:solid;border-width:1px;overflow:hidden;word-break:normal;}
% .tg .tg-rkwu{font-family:"Courier New", Courier, monospace !important;;vertical-align:top}
% .tg .tg-h31u{font-family:Arial, Helvetica, sans-serif !important;;vertical-align:top}
% .tg .tg-9hbo{font-weight:bold;vertical-align:top}
% .tg .tg-yw4l{vertical-align:top}
% </style>
% <table class="tg">
%   <tr>
%     <th class="tg-9hbo">Option</th>
%     <th class="tg-9hbo">Value</th>
%     <th class="tg-9hbo">Description</th>
%   </tr>
%   <tr>
%     <td class="tg-rkwu">'Scheme'</td>
%     <td class="tg-rkwu">'classic'</td>
%     <td class="tg-yw4l">Architecture or scheme of the Hopfield model</td>
%   </tr>
%   <tr>
%     <td class="tg-rkwu">'SimFcn'</td>
%     <td class="tg-rkwu">'euler' | 'runge-kutta' | 'talavan-yanez' (default)</td>
%     <td class="tg-yw4l">CHN Simulation method</td>
%   </tr>
%   <tr>
%     <td class="tg-rkwu">'TrainFcn'</td>
%     <td class="tg-rkwu">'traingty'</td>
%     <td class="tg-yw4l">Training method for the CHN applied to the GQKP</td>
%   </tr>
%   <tr>
%     <td class="tg-rkwu">'CheckpointPath'</td>
%     <td class="tg-rkwu">'' (default) | character vector</td>
%     <td class="tg-yw4l">Path to store the simulation trayectory</td>
%   </tr>
%   <tr>
%     <td class="tg-rkwu">'Dt'</td>
%     <td class="tg-rkwu">0.01 (default) | scalar positive value</td>
%     <td class="tg-yw4l">Integration step for 'euler' and 'runge-kutta' methods</td>
%   </tr>
%   <tr>
%     <td class="tg-rkwu">'E'</td>
%     <td class="tg-rkwu">13 (default) | scalar positive integer value</td>
%     <td class="tg-yw4l">Exponent for the stopping criteria tolerance, 1e-E</td>
%   </tr>
%   <tr>
%     <td class="tg-rkwu">'ExecutionEnvironment'</td>
%     <td class="tg-rkwu">'cpu' (default) | 'gpu'</td>
%     <td class="tg-yw4l">Execution environment of the CHN</td>
%   </tr>
%   <tr>
%     <td class="tg-rkwu">'MaxIter'</td>
%     <td class="tg-rkwu">2000 (default) | scalar positive integer value</td>
%     <td class="tg-yw4l">Maximum number of iterations</td>
%   </tr>
%   <tr>
%     <td class="tg-rkwu">'Q'</td>
%     <td class="tg-rkwu">0.8 (default) | scalar positive value</td>
%     <td class="tg-yw4l">Integration step reduction in 'talavan-yanez' during the first 'R_Iter' iterations</td>
%   </tr>
%   <tr>
%     <td class="tg-rkwu">'R_Iter'</td>
%     <td class="tg-rkwu">20 (default) | scalar positive integer value</td>
%     <td class="tg-yw4l">Number of iterations in 'talavan-yanez' where the integration step is reduced</td>
%   </tr>
%   <tr>
%     <td class="tg-rkwu">'SimulationPlot'</td>
%     <td class="tg-rkwu">false (default) | true</td>
%     <td class="tg-yw4l">Indicator to visualize the simulation</td>
%   </tr>
%   <tr>
%     <td class="tg-rkwu">'SimulationPlotPauseTime'</td>
%     <td class="tg-rkwu">0.8 (default) | scalar positive value</td>
%     <td class="tg-yw4l">Pause time in between iterations in the visualization of the simulation</td>
%   </tr>
%   <tr>
%     <td class="tg-rkwu">'TransferFcn'</td>
%     <td class="tg-rkwu">'satlin' | 'tanh' (default)</td>
%     <td class="tg-yw4l">Activation or Transfer function</td>
%   </tr>
%   <tr>
%     <td class="tg-rkwu">'U0'</td>
%     <td class="tg-rkwu">0.3 (default) | scalar positive value</td>
%     <td class="tg-yw4l">Slope of the transfer function</td>
%   </tr>
%   <tr>
%     <td class="tg-rkwu">'Verbose'</td>
%     <td class="tg-rkwu">false (default) | true</td>
%     <td class="tg-yw4l">Indicator to show information on the simulation of the CHN</td>
%   </tr>
%   <tr>
%     <td class="tg-rkwu">'VerboseFrequency'</td>
%     <td class="tg-rkwu">25 (default) | scalar positive integer value</td>
%     <td class="tg-yw4l">Frequency of verbose printing</td>
%   </tr>
% </table>
% </html>
%
%% 
% See <Example_GQKPusingCHN.html example> to solve a GQKP using a CHN. 