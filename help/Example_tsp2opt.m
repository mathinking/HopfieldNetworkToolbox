%% Improving Hopfield Networks performance using the CHN as a 2-opt
% In order to improve the performance of this heuristic technique,
% consecutive second phases of the Divide-and-Conquer scheme using 4 cities
% and 2 chains can be used, starting with an initial solution. This process 
% behaves like a 2-opt algorithm (in particular Lin-Kernighan's algorithm).

%% TSPLIB problem and network parameters
rng(6); % For reproducibility

%%
% TSPLIB Problem:
problem = tsplib({'berlin52'});

%%
% Number of cities:
N = problem.NumberOfCities;

%%
% Free parameter C:
C = 1e-5;

%% Creating the |HopfieldNetworkTSP| object using the Divide-and-Conquer simulation method
% Providing problem coordinates cities and distance matrix to the
% |HopfieldNetworkTSP| network by creating a |HopfieldNetworkTSPObject| 
% object of options
options = tsphopfieldnetOptions('Coordinates',problem.Coordinates,...
                                'DistanceMatrix',problem.DistanceMatrix,...
                                'DistanceType',problem.DistanceType,...
                                'Scheme','classic&2opt');
net = tsphopfieldnet(N,C,options);

%% Training the network
% The default training algorithm is |trainty|
train(net);

%% Simulating the network
% The simulation is using the algorithm is |talavan-yanez|
sim(net);

%% Visualizing results
getResults(net)
plot(net);
