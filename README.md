# Hopfield Network Toolbox for MATLAB
This Hopfield Network Toolbox is mainly focused in **Continuous Hopfield Networks** (CHNs). This Toolbox is based on the work by Javier Y��ez, Pedro M. Talav�n and Lucas Garc�a.

The Continuous Hopfield Network (CHN) is a **recurrent neural network** with an associated differential equation, whose state evolves from an initial condition to an equilibrium point by minimizing a Lyapunov function. 
As the Lyapunov function is associated with an objective function of the optimization problem (i.e. the mapping process), the equilibrium, or stable point, helps identify a local optimum for the optimization problem.
The dynamics of the CHN is described by a **differential equation**:

![alt tag](http://mathurl.com/hzhnzj5.png) 

and the output function is a hyperbolic tangent:

![alt tag](http://mathurl.com/zdeg52h.png)

The existence of an equilibrium point is guaranteed if a Lyapunov or energy function exists. The idea is that the network's Lyapunov function, when ![alt tag](http://mathurl.com/hl2t4by.png), is associated with the cost function to be minimized in the combinatorial problem.

The CHN will solve those combinatorial problems which can be expressed as the constrained minimization of: 

![alt tag](http://mathurl.com/goa77to.png)

However, at this point the Hopfield Network Toolbox is primarily designed to solve the Traveling Salesman Proble,.

## Download 
It is recommended to use the [latest release](https://github.com/mathinking/HopfieldNetworkToolbox/releases/latest). You may download the entire source code or a single installable Toolbox file.

## Install
To install the Toolbox, simply run

```sh
>> setup_hopfieldNetwork
```
in MATLAB's Command Window.

## Development
Want to contribute? Great! Feel free to fork the repository and contact us for instructions and suggestions. 

## Questions?
Open a new Issue and label it as a question. We will get back to you.

## References
  - [�Neural� computation of decisions in optimization problems](http://www.ams.org/mathscinet-getitem?mr=824597)
  - [A continuous Hopfield network equilibrium points algorithm](http://www.sciencedirect.com/science/article/pii/S0305054804000243)
  - [Parameter setting of the Hopfield network applied to TSP](http://www.sciencedirect.com/science/article/pii/S0893608002000217)
  - [Improving the Hopfield model performance when applied to the traveling salesman problem: A divide-and-conquer scheme](http://link.springer.com/article/10.1007/s00500-016-2039-8)

## Contact us
Send us an [Email](mailto:lucasgarciarodriguez@ucm.es) with your comments/suggestions.

## License
BSD 2-clause �Simplified� License
