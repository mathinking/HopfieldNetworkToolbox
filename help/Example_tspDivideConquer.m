%% Improving Hopfield Networks performance using a Divide-and-Conquer scheme
% In order to improve the performance of this heuristic technique, a 
% Divide-and-Conquer strategy based on two phases is proposed. The first
% phase involves linking cities with the most neighbors to define a set of
% chains of cities and, secondly, to join these with isolated cities to 
% define the final tour. Both problems are solved by mapping the two TSPs 
% onto their respective CHNs.

%% TSPLIB problem and network parameters
rng(3); % For reproducibility

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
                                'Scheme','divide-conquer',...
                                'Tau',2);
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
