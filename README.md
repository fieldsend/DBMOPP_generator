# DBMOPP_generator

Codebase for Distance-based Multi-Objective Point Problem instances.

This repository holds Matlab code for the Distance-Based Multi-Objective Point Problem (DBMOPP) instance generator. 

The DBMOPP generator in the <code>DBMOPP.m</code> file is the most recent version of the generator, and follows an object-oriented design. 

This is a substantial new version of the generator including the integration of hard and soft constraints and additional features, recently acceptanced in IEEE Transactions on Evolutionary Computation.

Jonathan E. Fieldsend, Tinkle Chugh, Richard Allmendinger, and Kaisa Miettinen. 
2021. A Visualizable Test Problem Generator for Many-Objective Optimization, 
IEEE Transactions on Evolutionary Computation, 
to appear. doi: https://10.1109/TEVC.2021.3084119.

All instance generation functionality is in this single file. Note there are a few "quality of life" methods still still be completed over the next month which will allow you to e.g. sample Pareto members, and some additionally plotting functions (the stubs for these are in the class and details are below). Illustrations of usage directly. Details regarding the earlier version of the generator can be found towards the bottom of this page.

## DBMOPP class

### Public methods
The public methods of this class are as follows:

<code>DBMOPP(numberOfObjectives, numberOfDesignVariables, numberOfLocalParetoSets, numberOfDominanceResistanceRegions, numberOfGlobalParetoSets, proportionOfConstrainedSpaceIfChecker, globalParetoSetType, constraintType, numberOfdiscontinousObjectiveFunctionRegions, variableSolutionDensity, varyingObjectiveScales, proportionOfNeutralSpace, monte_carlo_samples)</code> Which constructs an instance of the problem 

<code>plotLandscapeForSingleObjective(index, resolution)</code> which plots the 2D landscape of the objective numbered by <code>index</code>. Plot generated by a grid with <code>resolution</code> samples at each dimension.

<code>plotParetoSetMembers(resolution)</code> Samples points on a grid with <code>resolution</code> samples at each dimension, plots those which are Pareto optimal for this instance.

<code>plotDominanceLandscape(resolution)</code> Plots the dominance landscape of this problem, by generating by a grid with <code>resolution</code> samples at each dimension.

<code>plotProblemInstance()</code> Plots the different regions set up for this problem instance.

<code>isAParetoSetMember(x, suppressWarning)</code> Returns true if <code>x</code> is a Pareto set member, false otherwise. Obviously you should not be using this method during an optimisation! Pass <code>suppressWarning</code> as true if you don't want to be reminded of this...

<code>getAParetoSetMember()</code> Returns a Pareto decision vector uniformly at random for this instance.

<code>evaluate(x)</code> evaluates a design vector <code>x</code>

Help information can be accessed at the command line for each of the public methods, e.g. the command

<code>help DBMOPP/evaluate</code>

will display

```
    [objective_vector, soft_constraint_violation, hard_constraint_violation] = evaluate(obj,x)
 
    Evalutes the design vector x under this instance of the problem
 
    INPUTS
    x = design vector to evaluate, should adhere to the box
      constraints of the problem, i.e. between -1 and +1 on all
      dimensions
 
    OUTPUTS
  
    objective_vector = vector of objective values associated with
      x for this DBMOPP problem. If x violates any soft or hard
      constraints will be a vector of NaNs
    soft_constraint_violation = Soft constraint violation. 0 if
      no violation, otherwise value increases with distance to the
      boundary from the constrain/legal space
      hard_constraint_violation = Hard constraint violation. Value
      of true if there is a violation, false otherwise 
```

### Example usage

First let's create an instance from the generator, this is done with the constructor, which creates a DBMOPP instance based on the arguments. Default values are used when arguments are missing (see documentation in code). An example would be

```Matlab
my_instance = DBMOPP(4,2,0,0,5,0,1,0,0,false,false,0);
```
This creates <code>my_instance</code> which has 4 objectives, 2 decision variables, 5 disconnected Pareto set regions which have global Pareto set type '1', meaning they are partially intersecting -- the entire Pareto front can be described by fewer than five of the regions, but not one alone.

Calling 
```Matlab
my_instance.plotProblemInstance();
```
plots a helpful visualisation of the problem as constructed, in this case, 
 
![Constructed problem instance](/assets/images/instance_regions.jpg "Constructed problem instance")

As we don't have e.g. any constrained space, or neutral space in this instance example, it is showing the attractor locations (labelled with the objective they minimise) and the convex hull of the region they bound.
 
Calling <code>my_instance.plotParetoSetMembers()</code> plots which samples on the default resolution (a 500 by 500 grid) are Pareto optimal. As we have set up a partially intersecting Pareto set type, some areas in the convex hull are additionally penalised (the objective values increased), meaning they are not Pareto optimal, this plot shows us the Pareto optimal locations from the grid of samples 

![Pareto optimal points from grid](/assets/images/instance_pareto.jpg "Pareto optimal points from grid")
 
If we look at the fitness landscape for objectve 1, by calling 
```Matlab
my_instance.plotLandscapeForSingleObjective(1);
``` 
we can see how this is affected both by the attractors in each region, and the offset applied inside the convex hull of the attractor groupings defining the regions, to induce the partially overlapping Pareto sets. 
 
![Objective 1 landscape](/assets/images/example_objective1.jpg "Objective 1 landscape")

The <code>evaluate</code> method allows you to utilise the generated instance to assess the quality of a design, it returns the objective vector, and the soft and hard constrain violation values. E.g.
    
```Matlab
x = [0.6 0.4];
[y, sc, hc] = my_instance.evaluate(x);
```
Evaluates the design <code>x</code>. This instance has no constraints (bar box constraints), so the <code>sc</code> and <code>hc</code> values are both zero.

```Matlab
y
y =

    0.3218    0.1977    0.5040    0.4860
sc

sc =

     0

hc

hc =

     0
 
```    

The number of solutions found which are Pareto optimal under the instance configuration can be assessed using the <code>isAParetoSetMember</code> method:

```Matlab
my_instance.isAParetoSetMember(x,true)

ans =

  logical

   0
```
In this case, <code>x</code> is not Pareto optimal (as can also be seen in the earlier plots), a close neighbour is though:

```Matlab
my_instance.isAParetoSetMember([-0.487,0.4],true)

ans =

  logical

   1
```

Let us generate a more complicated instance

```Matlab
my_complex_instance = DBMOPP(3,10,3,3,2,0.3,2,8,2,false,true,0.2);
```
This is a 10-dimensional problem, with three competing objectives. 30% of the design space is soft constraint violating, 20% is neutral, there are 2 global Pateto set regions, which do not duplicate performance. There are also local fronts and domainace resistance regions.

![Complex problem instance](/assets/images/complex_dbmopp.jpg "Complex problem instance")

### (Current) stub methods

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

