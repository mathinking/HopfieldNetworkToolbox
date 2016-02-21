%% Hopfield Network Toolbox Release Notes
% This toolbox is solves combinatorial optimization problems by using
% Continuous Hopfield Networks (CHN). 
%
% The Hopfield Network Toolbox requires the following products: HERE
% 
%% Version 1.1.1, 21-Feb-2015
% Bug fixing for Unix platforms
%
% List of <https://github.com/mathinking/HopfieldNetworkToolbox/issues?q=milestone%3Av1.1.1+is%3Aclosed+label%3Abug
% Bugs> fixed in this release.
% 
% <https://github.com/mathinking/HopfieldNetworkToolbox/releases/tag/v1.1.1
% Version 1.1.1 on GitHub> 

%% Version 1.1, 21-Feb-2015
% New App, Toolbox Documentation and Examples.
% 
% <https://github.com/mathinking/HopfieldNetworkToolbox/issues?q=milestone%3Av1.1+is%3Aclosed+label%3Aenhancement 
% New Features>:
%%
% * Hopfield Net TSP solver App
% * Toolbox documentation
% * Step-by-step examples
% * TSPLIB automatic download
%
% List of <https://github.com/mathinking/HopfieldNetworkToolbox/issues?q=milestone%3Av1.1+is%3Aclosed+label%3Abug
% Bugs> fixed in this release.
% 
% <https://github.com/mathinking/HopfieldNetworkToolbox/releases/tag/v1.1
% Version 1.1 on GitHub> 

%% Version 1.0, 02-Nov-2015
% Initial release of Hopfield Network Toolbox.
% 
% This release is mainly focused in solving the Traveling Salesman Problem 
% using the Continuous Hopfield Network (_CHN_). However, the release also 
% provides a class structure to solve generic combinatorial optimization
% problems. Development in this area is undergoing.
% 
% The class to solve the TSP problems using _CHNs_ is |tsphopfieldnet|. This 
% network can solve any TSP problem, provided its coordinates or distance 
% matrix. 
% The Toolbox also includes the library TSPLIB, a de facto library for TSP
% benchmarks. Instances with up to 13509 cities have been tested using the  
% |tsphopfieldnet| network. Note that solving such instances might require 
% a large amount of memory.
% 
% Two main algorithms can be tested in this release:
% 
% * |talavan-yanez|: based on the paper 
%   <http://www.sciencedirect.com/science/article/pii/S0893608002000217 
%   Parameter setting of the Hopfield network> applied to TSP by _Pedro M. 
%   Talaván and Javier Yáñez_.
% * |divide-conquer|: based on the paper (pending publishing at the time of
%   this release) Improving the Hopfield model performance when applied to 
%   the traveling salesman problem: A divide-and-conquer scheme by _Lucas 
%   García, Pedro M. Talaván and Javier Yáñez_.
%
% <https://github.com/mathinking/HopfieldNetworkToolbox/releases/tag/v1.0
% Version 1.0 on GitHub> 
%
%% References
% * García, L., Talaván, P. M., & Yáñez, J. (2016). 
% <http://link.springer.com/article/10.1007/s00500-016-2039-8 Improving the
% Hopfield model performance when applied to the traveling salesman 
% problem>. Soft Computing, 1-15. 
% * Talaván, P. M., & Yáñez, J. (2002). 
% <http://www.sciencedirect.com/science/article/pii/S0893608002000217 
% Parameter setting of the Hopfield network applied to TSP>. Neural 
% Networks, 15(3), 363-373.
% * Talaván, P. M., & Yáñez, J. (2005). 
% <http://www.sciencedirect.com/science/article/pii/S0305054804000243 A 
% continuous Hopfield network equilibrium points algorithm>. Computers & 
% operations research, 32(8), 2179-2196.
% 
% The latest version of this toolbox is available on
% <https://github.com/mathinking/HopfieldNetworkToolbox/releases/latest 
% GitHub>. 
%
