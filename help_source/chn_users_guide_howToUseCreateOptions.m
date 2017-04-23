%% How to use tsphopfieldnet and createOptions
% |tsphopfieldnet| is the class that creates Continuous Hopfield Networks
% for solving the Traveling Salesman Problem.
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
% Optionally, |tsphopfieldnet| may receive a third input argument that
% gives further information of the TSP problem considered. This should be
% the case for most of the TSP problems you may want to solve.
%
%%
% *Syntax*
%
%      net = tsphopfieldnet(N,C)
%      net = tsphopfieldnet(N,C,options)
% 
%% How |tsphopfieldnetOptions| work
% |tsphopfieldnetOptions| is a class for |tsphopfieldnet| that allows to specify
% the details of each TSP problem.
% 
%%
% *Syntax*
%
%      options = tsphopfieldnetOptions('param1','value1','param2','value2',...)

%%
% The following table lists the available options for |tsphopfieldnet|.
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
%     <td class="tg-rkwu">'Coordinates'</td>
%     <td class="tg-rkwu">Nx2 real numeric matrix</td>
%     <td class="tg-yw4l">TSP coordinates</td>
%   </tr>
%   <tr>
%     <td class="tg-rkwu">'DistanceMatrix'</td>
%     <td class="tg-rkwu">NxN real numeric matrix</td>
%     <td class="tg-yw4l">TSP distance matrix</td>
%   </tr>
%   <tr>
%     <td class="tg-rkwu">'DistanceType'</td>
%     <td class="tg-rkwu">'EUC'<br>'EUC_2D'<br>'GEO'<br>'ATT'<br>'CEIL_2D'<br>'EXPLICIT'</td>
%     <td class="tg-yw4l">Distance type</td>
%   </tr>
%   <tr>
%     <td class="tg-rkwu">'Names'</td>
%     <td class="tg-rkwu">1xN cell</td>
%     <td class="tg-yw4l">TSP city names</td>
%   </tr>
%   <tr>
%     <td class="tg-rkwu">'Subtours'</td>
%     <td class="tg-rkwu">1xN cell</td>
%     <td class="tg-yw4l">TSP fixed subtours</td>
%   </tr>
%   <tr>
%     <td class="tg-rkwu">'SubtoursPositions'</td>
%     <td class="tg-rkwu">1xN matrix</td>
%     <td class="tg-yw4l">TSP fixed subtours starting positions</td>
%   </tr>
%   <tr>
%   <tr>
%     <td class="tg-rkwu">'Tau'</td>
%     <td class="tg-rkwu">1x1 double</td>
%     <td class="tg-yw4l">Number of neighboring cities</td>
%   </tr>
%     <td class="tg-rkwu">'SimFcn'</td>
%     <td class="tg-rkwu">'euler'<br>'talavan-yanez'<br>'divide-conquer'</td>
%     <td class="tg-h31u">Simulation method</td>
%   </tr>
%   <tr>
%     <td class="tg-rkwu">'TrainFcn'</td>
%     <td class="tg-rkwu">'trainty'</td>
%     <td class="tg-yw4l">Training method</td>
%   </tr>
%   <tr>
%     <td class="tg-rkwu">'K'</td>
%     <td class="tg-rkwu">1x1 double</td>
%     <td class="tg-yw4l">Number of chains</td>
%   </tr>
%   <tr>
%     <td class="tg-rkwu">'CheckpointPath'</td>
%     <td class="tg-rkwu">char</td>
%     <td class="tg-yw4l">Valid writable path to store simulation</td>
%   </tr>
%   <tr>
%     <td class="tg-rkwu">'Dt'</td>
%     <td class="tg-rkwu">1x1 double</td>
%     <td class="tg-yw4l">Step for Euler method</td>
%   </tr>
%   <tr>
%     <td class="tg-rkwu">'E'</td>
%     <td class="tg-rkwu">1x1 double</td>
%     <td class="tg-yw4l">Tolerance exponent</td>
%   </tr>
%   <tr>
%     <td class="tg-rkwu">'ExecutionEnvironment'</td>
%     <td class="tg-rkwu">'CPU'<br>'GPU'</td>
%     <td class="tg-yw4l">Hardware resources</td>
%   </tr>
%   <tr>
%     <td class="tg-rkwu">'MaxIter'</td>
%     <td class="tg-rkwu">1x1 double</td>
%     <td class="tg-yw4l">Maximum number of iterations</td>
%   </tr>
%   <tr>
%     <td class="tg-rkwu">'Q'</td>
%     <td class="tg-rkwu">1x1 double</td>
%     <td class="tg-yw4l">Parameter factor for reducing <br>integration step</td>
%   </tr>
%   <tr>
%     <td class="tg-rkwu">'R_Iter'</td>
%     <td class="tg-rkwu">1x1 double</td>
%     <td class="tg-yw4l">Number of iterations for reducing<br>integration step</td>
%   </tr>
%   <tr>
%     <td class="tg-rkwu">'SimulationPlot'</td>
%     <td class="tg-rkwu">1x1 logical</td>
%     <td class="tg-yw4l">Show simulation process plot</td>
%   </tr>
%   <tr>
%     <td class="tg-rkwu">'SimulationPlotPauseTime'</td>
%     <td class="tg-rkwu">1x1 double</td>
%     <td class="tg-yw4l">Pause time during simulation <br> process plot</td>
%   </tr>
%   <tr>
%     <td class="tg-rkwu">'TransferFcn'</td>
%     <td class="tg-rkwu">'satlin'<br>'tanh'</td>
%     <td class="tg-yw4l">Transfer function</td>
%   </tr>
%   <tr>
%     <td class="tg-rkwu">'U0'</td>
%     <td class="tg-rkwu">1x1 double</td>
%     <td class="tg-yw4l">Initial potential</td>
%   </tr>
%   <tr>
%     <td class="tg-rkwu">'Verbose'</td>
%     <td class="tg-rkwu">1x1 logical</td>
%     <td class="tg-yw4l">Display information of <br>simulation process</td>
%   </tr>
%   <tr>
%     <td class="tg-rkwu">'VerboseFrequency'</td>
%     <td class="tg-rkwu">1x1 integer value</td>
%     <td class="tg-yw4l">Frequency of verbose printing</td>
%   </tr>
% </table>
% </html>
%
%% 
% See Examples using <Example_tspUsingTSPLIB.html TSPLIB> or when 
% <Example_tspUsingCoords.html city coordinates> or 
% <Example_tspUsingDistance.html distance matrix> are provided.
