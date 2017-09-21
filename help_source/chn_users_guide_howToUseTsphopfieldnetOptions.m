%% How to use tsphopfieldnet and tsphopfieldnetOptions
% |tsphopfieldnet| is a function that returns an object of the class
% |HopfieldNetworkTSP|. This object is a Continuous Hopfield Network
% designed to solve the Traveling Salesman Problem.
% 
% At the same time, |tsphopfieldnetOptions| returns an object of the class
% |HopfieldNetworkTSPOptions|. Such object may be used to customize the
% Hopfield network object from |HopfieldNetworkTSP|.
%
%% How |tsphopfieldnet| works
% |tsphopfieldnet|'s default behavior is to take only 2 input arguments:
%%
% * N: number of cities
% * C: free parameter (greater than 0)
% 
% In this case, the cities are place in the vertices of regular polygons
% (see <Example_tspUsingRegularPolygons.html Example>).
% 
% Optionally, |tsphopfieldnet| may receive a third input argument (of the 
% class |HopfieldNetworkTSPOptions| using the function 
% |tsphopfieldnetOptions|) which gives further information of the TSP 
% problem considered. This way an arbitrary TSP problem can be created. 
%
%%
% *Syntax*
%
%      net = tsphopfieldnet(N,C)
%      net = tsphopfieldnet(N,C,options)
% 
%% How |tsphopfieldnetOptions| works
% |tsphopfieldnetOptions| is a function returning an object of the class 
% |HopfieldNetworkTSPOptions| which allows to specify all the settings and 
% details of the TSP problem.
% 
%%
% *Syntax*
%
%      options = tsphopfieldnetOptions('param1','value1','param2','value2',...)

%%
% The following table lists the available options for 
% |HopfieldNetworkTSPOptions|.
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
%     <td class="tg-rkwu">'classic' (default) | 'divide-conquer' | '2opt' | 'classic&2opt' | 'divide-conquer&2opt'</td>
%     <td class="tg-yw4l">Architecture or scheme of the Hopfield model</td>
%   </tr>
%   <tr>
%     <td class="tg-rkwu">'SimFcn'</td>
%     <td class="tg-rkwu">'euler' | 'runge-kutta' | 'talavan-yanez' (default)</td>
%     <td class="tg-yw4l">CHN Simulation method</td>
%   </tr>
%   <tr>
%     <td class="tg-rkwu">'TrainFcn'</td>
%     <td class="tg-rkwu">'trainty'</td>
%     <td class="tg-yw4l">Training method for the CHN applied to the TSP</td>
%   </tr>
%   <tr>
%     <td class="tg-rkwu">'K'</td>
%     <td class="tg-rkwu">0 (default) | scalar integer value greater or equal than 0</td>
%     <td class="tg-yw4l">Number of fixed chains</td>
%   </tr>
%   <tr>
%     <td class="tg-rkwu">'Coordinates'</td>
%     <td class="tg-rkwu">[] (default) | real-value 2-column matrix</td>
%     <td class="tg-yw4l">TSP Cities' coordinates</td>
%   </tr>
%   <tr>
%     <td class="tg-rkwu">'DistanceMatrix'</td>
%     <td class="tg-rkwu">[] (default) | real-value square matrix</td>
%     <td class="tg-yw4l">TSP Distance Matrix</td>
%   </tr>
%   <tr>
%     <td class="tg-rkwu">'DistanceType'</td>
%     <td class="tg-rkwu">'geo' | 'euc_2d' | 'euc' (default) | 'att' | 'ceil_2d' | 'explicit'</td>
%     <td class="tg-yw4l">Distance type between TSP cities</td>
%   </tr>
%   <tr>
%     <td class="tg-rkwu">'Names'</td>
%     <td class="tg-rkwu">'' (default) | cell array of strings</td>
%     <td class="tg-yw4l">TSP city names</td>
%   </tr>
%   <tr>
%     <td class="tg-rkwu">'PlotPhase'</td>
%     <td class="tg-rkwu">false (default) | true</td>
%     <td class="tg-yw4l">Indicator to display the result of the phases in 'divide-conquer'</td>
%   </tr>
%   <tr>
%     <td class="tg-rkwu">'Subtours'</td>
%     <td class="tg-rkwu">'' (default) | cell array of strings</td>
%     <td class="tg-yw4l">Chains of cities separated by hyphens</td>
%   </tr>
%   <tr>
%     <td class="tg-rkwu">'SubtoursPositions'</td>
%     <td class="tg-rkwu">[] (default) | vector of integer values with as many columns as elements in 'Subtours'</td>
%     <td class="tg-yw4l">Location of the first city in each pre-fixed chains of cities in 'Subtours'</td>
%   </tr>
%   <tr>
%     <td class="tg-rkwu">'Tau'</td>
%     <td class="tg-rkwu">[] (default) | scalar integer value</td>
%     <td class="tg-yw4l">Number of neighboring cities considered in the first phase of 'divide-conquer'</td>
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
% See Examples using <Example_tspUsingTSPLIB.html TSPLIB> or when 
% <Example_tspUsingCoords.html city coordinates> or 
% <Example_tspUsingDistance.html distance matrix> are provided.
