%% Solving the TSP with TSPLIB cities using a Hopfield Network
% TSPLIB is an extensive set of problems frequently used for solving the 
% TSP (and related problems) and used as a powerful reference for
% benchmarking different algorithms. 
% In this section, |tsphopfieldnet| is used for solving the TSP when the 
% problem coordinates are given by one of the library problems from TSPLIB. 
% 

%% TSPLIB problem and network parameters
rng(6); % For reproducibility

%%
% TSPLIB Problem:
problem = tsplib({'berlin52'});

%%
% Number of cities:
N = problem.nCities;

%%
% Free parameter C:
C = 0.00001;
 
%% Creating the |tsphopfieldnetwork| object
% Providing problem coordinates cities and distance matrix to the
% |tsphopfieldnet| network by creating a structure of options
options = tsphopfieldnet.createOptions('coords',problem.coords,'d',problem.d,'type',problem.type);
net = tsphopfieldnet(N,C,options);

%%
% Data coordinates (cities) can be visualized before training takes place:
plot(net);

%% Training the network
% The default training algorithm is |trainty|
train(net);
%%
% Results of the training phase. Network parameters
getTrainParam(net)

%% Simulating the network
% The default simulation algorithm is |talavan-yanez|
sim(net);

%% Visualizing results
getResults(net)
plot(net);
