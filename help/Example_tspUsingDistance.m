%% Solving the TSP provided the distance matrix using Hopfield Network
% The following example shows how to solve any TSP problem, provided the
% distance matrix between cities.
%
%% Network parameters
rng(22); % For reproducibility

%%
% Free parameter C:
C = 0.00001;
 
%% Creating the |tsphopfieldnetwork| object
% Providing the distance matrix
d = [.0000 .3361 .3141 .3601 .5111 .5176 .2982 .4564 .3289 .2842;...
     .3361 .0000 .1107 .6149 .8407 .8083 .5815 .6418 .4378 .3934;...
     .3141 .1107 .0000 .5349 .7919 .8207 .5941 .6908 .4982 .4501;...
     .3601 .6149 .5349 .0000 .3397 .6528 .5171 .7375 .6710 .6323;...
     .5111 .8407 .7919 .3397 .0000 .4579 .4529 .6686 .7042 .6857;...
     .5176 .8083 .8207 .6528 .4579 .0000 .2274 .2937 .4494 .4654;...
     .2982 .5815 .5941 .5171 .4529 .2274 .0000 .2277 .2690 .2674;...
     .4564 .6418 .6908 .7375 .6686 .2937 .2277 .0000 .2100 .2492;...
     .3289 .4378 .4982 .6710 .7042 .4494 .2690 .2100 .0000 .0498;...
     .2842 .3934 .4501 .6323 .6857 .4654 .2674 .2492 .0498 .0000];
 
%%      
% Providing the distance matrix to the network through the structure of 
% options      
options = tsphopfieldnet.createOptions('d',d);

%%
% Number of cities:
N = size(d,1);

net = tsphopfieldnet(N,C,options);

%% Training the network
% The default training algorithm is |trainty|
train(net);
%%
% Results of the training phase. Network parameters
getTrainParam(net)

%% Simulating the network
% The default simulation algorithm is |talavan-yanez|
sim(net);

%% Tour Length
getResults(net,'tourLength')
