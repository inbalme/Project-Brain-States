
This folder contains Matlab and C code to generate figures from the manuscript:

The spatial structure of correlated neuronal variability.  R. Rosenbaum, M.A. Smith, A. Kohn, J.E. Rubin and B. Doiron. Nature Neuroscience, 2016.

Please cite this paper if you use any of the included software or any modifications of it.

The GPFA code in the subfolder gpfa-master was downloaded directly from: http://toliaslab.org/publications/ecker-et-al-2014/
and comes from the algorithms describes in the following papers:

Yu, B, et al. "Gaussian Process Factor Analysis for Low-Dimensional Single-Trial Analysis of Neural Population Activity." Journal of Neurophysiology 102.3 (2009).

Ecker, Alexander S., et al. "State dependence of noise correlations in macaque primary visual cortex." Neuron 82.1 (2014): 235-248.

Please cite these papers if you use the GPFA algorithms.

------

Instructions for use:
- Before running any of the Matlab code, first you must compile all of the C code using the mex compiler in Matlab.  Specifically, in the Matlab command line, run the following commands:
mex EIF2DSpatialNetwork.c
mex CalcRheoBaseEIF.c
mex EIFTwoPopNetwork.c
mex EIF2DSpatialNetworkGainMod.c

- After compiling the C code, you can run the network simulations by running the scripts NetworkSimForFigureX.m where X is the Figure number.
- After running each simulation, you can generate the figures by running the scripts MakeFigureX.m.  You need to run the network simulation once before you can generate the associated figure.


- Apologies for the sparse comments and occasionally confusing algorithms.  For questions about the code, please e-mail Robert Rosenbaum (currently at Robert.Rosenbaum@nd.edu)
