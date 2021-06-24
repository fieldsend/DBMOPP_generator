# DBMOPP_generator

Codebase for Distance-based Multi-Objective Point Problem instances.

This repository holds Matlab code for the Distance-Based Multi-Objective Point Problem (DBMOPP) instance genetor. 

The DBMOPP generator in the DBMOPP.m file is the most recent version of the generator, and follows an object-oriented design. 

This substaional new version of the generator, including the integration of hard and soft constraints and additional features, recently acceptanced in IEEE Transactions on Evolutionary Computation

Jonathan E. Fieldsend, Tinkle Chugh, Richard Allmendinger, and Kaisa Miettinen. 
2021. A Visualizable Test Problem Generator for Many-Objective Optimization, 
IEEE Transactions on Evolutionary Computation, 
to appear. doi: https://10.1109/TEVC.2021.3084119.

All instance generation functionality is in this version, there are a few "quality of life" methods still still be completed over the next month which allow you to sample Pareto members, and some additionally plotting functions, but the stubs are in the class and details are below, with illustrations of usage.

## DBMOPP class

First let's create an instance from the generator, this is done with the constructor:

 <code>DBMOPP(numberOfObjectives, numberOfDesignVariables, numberOfLocalParetoSets, numberOfDominanceResistanceRegions, numberOfGlobalParetoSets, proportionOfConstrainedSpaceIfChecker, globalParetoSetType, constraintType, numberOfdiscontinousObjectiveFunctionRegions, variableSolutionDensity, varyingObjectiveScales, proportionOfNeutralSpace, monte_carlo_samples)</code> 
 
creates a DBMOPP instance based on the arguments. Default values are used when arguments are missing (see documentation in code). An example would be

<code>my_instance = DBMOPP(4,2,0,0,5,0,1,0,0,false,false,0)</code>

This creates <code>my_instance</code> which has 4 objectives, 2 descision variables, 5 disconnected Pareto set regions which have global Pareto set type '1', meaning they are partially intersecting -- the entire Pareto front can be described by fewer than five of the regions, but not one alone.

Calling <code>my_instance.plotProblemInstance()</code> plots a helpful visualisation of the problem as constructed, in this case, 
 
![Constructed problem instance](/assets/images/instance_regions.jpg "Constructed problem instance")

As we don't have e.g. any constrained space, or neutral space in this instance example, it is showing the attractor locations (labelled with the objective they minimise) and the convex hull of the region they bound.
 
Calling <code>my_instance.plotParetoSetMembers()</code> plots which samples on the default resolution (a 500 by 500 grid) are Pareto optimal. As we have set up a partially intersecting Pareto set type, some areas in the convex hull are additionally penalised (the objective values increased), meaning they are not Pareto optimal, this plot shows us the Pareto optimal locations from the grid of samples 

![Pareto optimal points from grid](/assets/images/instance_pareto.jpg "Pareto optimal points from grid")
 
### Stub methods

The current methods are stubs and will throw errors, they will be filled in shortly

<code>plotDominanceLandscape</code> -- will plot the dominance landscape of an instance (exists in the older generator, just need to port for the OO version)

<code>getAParetoSetMember</code> -- will return, uniformly at random, a Pareto set member

## Older GECCO 2019 generator

Previous version of the generator with a subset of the functionality is described here:

Jonathan E. Fieldsend, Tinkle Chugh, Richard Allmendinger, and Kaisa Miettinen.
2019. A Feature Rich Distance-Based Many-Objective Visualisable
Test Problem Generator. In Genetic and Evolutionary Computation Conference
(GECCO ’19), July 13–17, 2019, Prague, Czech Republic. ACM, New York

Paper availabale from

Insititutional repository: https://ore.exeter.ac.uk/repository/handle/10871/36824

Publisher Doi: https://doi.org/10.1145/3321707.3321727

The codebase for this earlier paper is in the GECCO2019 folder.

USAGE:

Two main functions are 

distance_problem_generator -- this generates a problem instance (please see its help documentation)

distance_points_problem -- this takes a design vector and a problem instance and returns its quality vector (please see its help documentation)

Other .m files are helper/plotter functions.

It is recommended currently that you download and use the Release 1.0.0 (as new features are being added currently so current version should be considered unstable).
