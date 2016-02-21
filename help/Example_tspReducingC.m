%% Improving the convergence of the Hopfield Network applied to the TSP by reducing the parameter C
% In this section, you will be able to see how the convergence of the
% Hopfield Network is improved when reducing the free parameter C. Note
% that this had been proved experimentaly but there was no analytical proof
% until [ref].

%% TSPLIB problem and network parameters
seed = 10;
rng(seed); % For reproducibility

%%
% TSPLIB Problem:
problem = tsplib({'berlin52'});

%%
% Number of cities:
N = problem.nCities;

%%
% Free parameter C:
C = 1e5;
 
%% Creating the |tsphopfieldnetwork| object with value of C 'big'
% Providing problem coordinates cities and distance matrix to the
% |tsphopfieldnet| network by creating a structure of options
options = tsphopfieldnet.createOptions('coords',problem.coords,'d',problem.d,'type',problem.type);
net1 = tsphopfieldnet(N,C,options);

%%
% Data coordinates (cities) can be visualized before training takes place:
plot(net1);

%% Training the network
% The default training algorithm is |trainty|
train(net1);
%%
% Results of the training phase. Network parameters
getTrainParam(net1)

%% Simulating the network
% The default simulation algorithm is |talavan-yanez|
sim(net1);

%% Visualizing results
getResults(net1)
plot(net1);

%% Repeating the |tsphopfieldnet| network simulation with a smaller value of C

rng(seed); % For reproducibility
C = 1e-5;

net2 = tsphopfieldnet(N,C,options);
train(net2);
sim(net2);

getResults(net2)

plot(net2);

%%
% Note that not only the solution was improved, but the number of
% iterations and time to convergence was considerably improved.
