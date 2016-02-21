%% Concepts
%
%% What is a Neural Network?
%
% Artificial Neural Networks (ANN), or commonly Neural Networks, are a
% family of Machine Learning models, inspired by biological neural
% networks. ANNs were defined by Dr. Robert Hecht-Nielsen as:
%
% _"A computing system made up of a number of simple, highly interconnected   
% processing elements, which process information by their dynamic state
% response to external inputs."_
%
% You can find more information about Neural Networks online. One great
% introductory book is Neural Network Design, by _Martin T. Hagan, 
% Howard B. Demuth and Mark H. Beale_ (original authors of Neural Network 
% Toolbox). An ebook version is available from Martin T. Hagan's 
% <http://hagan.okstate.edu/nnd.html site>.
%
%% What is a Hopfield Neural Network?
%
% The Continuous Hopfield Network (CHN) is a *recurrent neural network* 
% with an associated differential equation, whose state evolves from an 
% initial condition to an equilibrium point by minimizing a Lyapunov 
% function. 
% As the Lyapunov function is associated with an objective function of the 
% optimization problem (i.e. the mapping process), the equilibrium, or 
% stable point, helps identify a local optimum for the optimization 
% problem.
%
% The dynamics of the CHN is described by a *differential equation*:
%%
% 
% $$\frac{du}{dt} = - \frac{u}{\Lambda} + T v + i^b$$
% 
% and the output function $v_i = g(u_i)$ is a hyperbolic tangent:
%%
% 
% $$g(u_i) = \frac{1}{2} \left( 1 + \tanh \left( \frac{u_i}{u_0} \right) \right), \qquad u_0 > 0$$
% 
%%
% The existence of an equilibrium point ($u^e$ such that 
% $u(t)=u^e \ \forall t \geq t_e$ for some $t_e \geq 0$) is guaranteed if a
% Lyapunov or energy function exists. The idea is that the network's 
% Lyapunov function, when $\Lambda \rightarrow{} \infty$, is associated 
% with the cost function to be minimized in the combinatorial problem.
%
% The CHN will solve those combinatorial problems which can be expressed as
% the constrained minimization of:
%%
%
% $$E(v) = -  \frac{1}{2} v^t T v - (i^b)^t v$$
%
%
%% The Traveling Salesman Problem
% Let $N$ be the number of cities in the TSP, and let $d_{xy}$ be the 
% distance between cities $x,y \in \{ 1,2,\ldots, N \}$.  Next, let $V$ be
% the $N \times N$ matrix of the state variable:
%
%%
% 
% $$v_{x,i} = \{ 1 \, \mbox{if the city} \, x \, \mbox{is visited in the 
% order} \, i \mbox{,} \, 0 \, \mbox{otherwise} \} $$
% 
% $V$ identifies a valid tour for the $TSP$ if the following constraints 
% are satisfied:
%
%%
% 
% * Every city must be visited only once: $S_x = \sum_{i=1}^N v_{x,i} = 1, 
% \quad \forall x \in \{1,2,\ldots,N \}$
% * Every position is associated with a unique city: $S_i = \sum_{x=1}^N 
% v_{x,i} = 1, \quad \forall i \in \{ 1,2,\ldots,N \}$
% 
% The *objective function* is: 
%
% $$ \min  \{ \frac{1}{2} \sum_{x=1}^N \sum_{y \neq x} \sum_{i=1}^N
% d_{xy} v_{xi}(v_{y(i+1)} + v_{y(i-1)}) \}, \mbox{(the $i+1$ 
% and $i-1$ subscripts are given modulo $N$)} $$ 
%
%% The Mapping Process
% 
% Given the state variable is $V\in [0,1]^{N\times N}$ for the TSP, the 
% energy function of the CHN is:
%%
%
% $$ E(V) = \frac{A}{2} \sum_x^N \sum_i^{N} \sum_{j \neq i}^{N} v_{x,i} 
% v_{x,j} + \frac{B}{2} \sum_i^{N} \sum_x^N \sum_{y \neq x}^N v_{x,i} 
% v_{y,i} + \frac{C}{2} ( \sum_x^N \sum_i^N v_{x,i} - N )^2 + \frac{D}{2} 
% \sum_x^N \sum_{y \neq x}^N \sum_i^N d_{x,y} v_{x,i} (v_{y,i-1} + 
% v_{y,i+1}) $$
% 
% By developing the four terms in this energy function, and comparing them 
% with the general form of the energy function, the weight of the arc 
% linking the neuron $(y,j)$ to $(x,i)$ and the incoming bias to neuron 
% $(x,i)$ are:
%%
%
% $$ T_{xi,yj} = -( A \delta_{x,y} (1-\delta_{i,j}) + B (1 - \delta_{x,y}) 
% \delta_{i,j} + C + D  d_{x,y} (\delta_{i,j-1} + \delta_{i,j+1}) ) $$
%
% $$i^b_{x,i} = C N $$
%
% with $x,y \in \{1,2,\ldots, N \}$ and $i,j \in \{1,2,\ldots, N \}$.
%
%% The Hopfield Network applied to the TSP
%%
% The Hopfield Network detailed in the previous section can be explained by
% using the following graph:
%
% <<network.png>>
% 
