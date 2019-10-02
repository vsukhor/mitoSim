#  MitoSim 

## A simulator of mitochondria reticulum (minimal spaceless formulation)

### Introduction

Mitochondria are organelles in the biological cells essential for keeping the cells alive. 
They are able to move around the cytosol (intracellular volume outside nucleus), where the form extensive quichly reorganizing networks.
Quantitative description oftheir morphology was problematic because of the seemingly irregular shape they can adopt. 

In the simulator, this problem is tackled by representing the mitochondria as a graph evolving in time. 
The graph dynamics results from the division(fission) and fusion reactions. 
As a result of the free energy minimization, in the coarse of the simulation run, the network settles in a dynamic steady-state configuration 
defined by the reaction parameters and the network dimensions.
A well-mixed spaceless representation is assumed, which means that any graph node has equal probability to interact with any other 
node allowed by the reaction.

For more detailes, please see the original publication:
Sukhorukov VM, Dikov D, Reichert AS, Meyer-Hermann M: Emergence of the mitochondrial reticulum from fission and fusion dynamics. 
PLoS Comput Biol. 2012, 8: e1002745-10.1371/journal.pcbi.1002745 [link]().

Please reference the above paper whenever results of this code are published.
