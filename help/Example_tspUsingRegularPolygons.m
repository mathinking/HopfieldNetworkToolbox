%% Solving the TSP with cities in Polygon Vertices using a Hopfield Network
% There are numerous TSP problems that have been used in the literature. In
% this section, an additional set of problem examples is proposed, by 
% placing TSP cities in the vertices of regular polygons. 
% This provides a whole new set of TSP problems that can be used to test
% different method approaches. 
%
% As a matter of fact, the coordinates of a regular polygon of size N is 
% used as the default behaviour of |tsphopfieldnet| if no set of 
% coordinates is provided to the network through |tsphopfieldnetOptions|.

%% Network parameters
rng(2); % For reproducibility

%%
% Number of cities:
N = 10;

%%
% Free parameter C:
C = 0.00001;
 
%% Creating the |HopfieldNetworkTSP| object
net = tsphopfieldnet(N,C);

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

%%
% The obtained tour is:
city(net, getResults(net,'VisitOrder'))
plot(net);
