%% Solving Continuous Hopfield Networks using Simulink

%% Opening the model
% To open the Continuous Hopfield Model in Simulink simply type in MATLAB:

open_system('continuousHopfieldNetwork')

%%
% The following example will solve the same GQKP proposed in a 
% <Example_GQKPusingCHN.html previous example>.

%% Setting-up the problem
P = [4,0;0,-2];
q = [0;0];
R = [1,1];
b = 1;

%% Building the weight matrix and bias vector from the obtained parametrization

%%
% Parametrization
alpha = 1;
Phi = 3;
eps = 3*alpha/2 + Phi/2;
beta = -alpha/2 - Phi;
Gamma(1,1) = (2*alpha + Phi/2);
Gamma(2,2) = Phi/2 - alpha;

T = -(alpha * P + R'*Phi*R - 2*Gamma);
ib = -(alpha * q + R'*beta + diag(Gamma));

%% Simulating the network
v0 = [0.6;0.2];
lambda = 1e5;
dt = 0.001;
output = sim('continuousHopfieldNetwork','StopTime','10');
output.yout(end,:)'

%%
% The obtained solution is in fact the global optimum.
