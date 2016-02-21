%% Improving the convergence of the Hopfield Network applied to the TSP by choosing a starting point close to the Saddle point
% In this section, you will be able to see how the convergence of the
% Hopfield Network is improved by choosing a starting point close to the 
% Saddle point of the Energy Function of the associated problem.
% 
%% Simulating the network with default starting point
seed = 5;
rng(seed); % For reproducibility

%%
% TSPLIB Problem:
problem = tsplib({'tsp225'});

%%
% Number of cities:
N = problem.nCities;

%%
% Free parameter C:
C = 1e-5;
 
%% 
% Creating the |tsphopfieldnetwork| object 
options = tsphopfieldnet.createOptions('coords',problem.coords,'d',problem.d,'type',problem.type);
net1 = tsphopfieldnet(N,C,options);

%% 
% Training the network
train(net1);

%% 
% Simulating the network
sim(net1);

%% Simulating the network with starting point close to Saddle point

net2 = tsphopfieldnet(N,C,options);
V = saddle(net2) + (rand(N) - 0.5)*1e-10;
invTransferFcn = getSetting(net2,'invTransferFcn');
U = invTransferFcn(V);
sim(net2,V,U);

%% Comparing results

getResults(net1,'tourLength')

getResults(net2,'tourLength')