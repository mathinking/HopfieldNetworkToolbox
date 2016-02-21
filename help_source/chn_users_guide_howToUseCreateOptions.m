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
%% How |createOptions| works
% |createOptions| is a method in |tsphopfieldnet| that allows to specify
% the details of each TSP problem.
% 
%%
% *Syntax*
%
%      options = tsphopfieldnet.createOptions('param1','value1','param2','value2',...)

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
%     <td class="tg-rkwu">'simFcn'</td>
%     <td class="tg-rkwu">'euler'<br>'talavan-yanez'<br>'divide-conquer'</td>
%     <td class="tg-h31u">Simulation method</td>
%   </tr>
%   <tr>
%     <td class="tg-rkwu">'trainFcn'</td>
%     <td class="tg-rkwu">'trainty'</td>
%     <td class="tg-yw4l">Training method</td>
%   </tr>
%   <tr>
%     <td class="tg-rkwu">'coords'</td>
%     <td class="tg-rkwu">Nx2 double</td>
%     <td class="tg-yw4l">TSP coordinates</td>
%   </tr>
%   <tr>
%     <td class="tg-rkwu">'d'</td>
%     <td class="tg-rkwu">NxN double</td>
%     <td class="tg-yw4l">TSP distance matrix</td>
%   </tr>
%   <tr>
%     <td class="tg-rkwu">'names'</td>
%     <td class="tg-rkwu">Nx1 cell</td>
%     <td class="tg-yw4l">TSP city names</td>
%   </tr>
%   <tr>
%     <td class="tg-rkwu">'type'</td>
%     <td class="tg-rkwu">'EUC'<br>'EUC_2D'<br>'GEO'<br>'ATT'<br>'CEIL_2D'<br>'EXPLICIT'</td>
%     <td class="tg-yw4l">Distance type</td>
%   </tr>
%   <tr>
%     <td class="tg-rkwu">'K'</td>
%     <td class="tg-rkwu">1x1 double</td>
%     <td class="tg-yw4l">Number of chains</td>
%   </tr>
%   <tr>
%     <td class="tg-rkwu">'tau'</td>
%     <td class="tg-rkwu">1x1 double</td>
%     <td class="tg-yw4l">Number of neighboring cities</td>
%   </tr>
%   <tr>
%     <td class="tg-rkwu">'hwResources'</td>
%     <td class="tg-rkwu">'CPU'<br>'GPU'</td>
%     <td class="tg-yw4l">Hardware resources</td>
%   </tr>
%   <tr>
%     <td class="tg-rkwu">'transferFcn'</td>
%     <td class="tg-rkwu">'satlin'<br>'tanh'</td>
%     <td class="tg-yw4l">Transfer function</td>
%   </tr>
%   <tr>
%     <td class="tg-rkwu">'invTransferFcn'</td>
%     <td class="tg-rkwu">'invsatlin'<br>'atanh'</td>
%     <td class="tg-yw4l">Inverse of transfer function</td>
%   </tr>
%   <tr>
%     <td class="tg-rkwu">'maxIter'</td>
%     <td class="tg-rkwu">1x1 double</td>
%     <td class="tg-yw4l">Maximum number of iterations</td>
%   </tr>
%   <tr>
%     <td class="tg-rkwu">'e'</td>
%     <td class="tg-rkwu">1x1 double</td>
%     <td class="tg-yw4l">Tolerance exponent</td>
%   </tr>
%   <tr>
%     <td class="tg-rkwu">'u0'</td>
%     <td class="tg-rkwu">1x1 double</td>
%     <td class="tg-yw4l">Initial potential</td>
%   </tr>
%   <tr>
%     <td class="tg-rkwu">'dt'</td>
%     <td class="tg-rkwu">1x1 double</td>
%     <td class="tg-yw4l">Step for Euler method</td>
%   </tr>
%   <tr>
%     <td class="tg-rkwu">'q'</td>
%     <td class="tg-rkwu">1x1 double</td>
%     <td class="tg-yw4l">Parameter for reducing <br>integration step</td>
%   </tr>
%   <tr>
%     <td class="tg-rkwu">'R_ITER'</td>
%     <td class="tg-rkwu">1x1 double</td>
%     <td class="tg-yw4l">Parameter for reducing<br>integration step</td>
%   </tr>
% </table>
% </html>
%
%% 
% See Examples using <Example_tspUsingTSPLIB.html TSPLIB> or when 
% <Example_tspUsingCoords.html city coordinates> or 
% <Example_tspUsingDistance.html distance matrix> are provided.
