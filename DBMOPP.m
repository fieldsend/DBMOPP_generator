% Class defining DBMOPP problem instance
%
% See: JE Fieldsend, T Chugh, R Allmendinger, K Miettinen (2021).
% A Visualizable Test Problem Generator for Many-Objective Optimization,
% IEEE Transactions on Evolutionary Computation,
% to appear. doi: 10.1109/TEVC.2021.3084119.
%
% Open source license details on github page for codebase:
%
% https://github.com/fieldsend/DBMOPP_generator
%
% Jonathan Fieldsend, University of Exeter, 2020, 2021

classdef DBMOPP < handle
    properties
        % arguments in
        numberOfObjectives % number of objectives >=2
        numberOfDesignVariables % number of design variables, >=2
        numberOfLocalParetoSets % number of local Pareto sets, >=0
        numberOfDominanceResistanceRegions % number of dominance resistance regions, >=0
        numberOfGlobalParetoSets % number of global Pareto sets, >= 1
        proportionOfConstrainedSpaceIfChecker % proprtion of space in 2D defined as constraint violating if checker type used 0 <= p < 1
        globalParetoSetType % 0, 1, 2 -- see help
        constraintType % 1, 2, 3, 4 -- see help
        numberOfdiscontinousObjectiveFunctionRegions % number of regions generated indcuing discontinous objective values on their border
        variableSolutionDensity % is solution density variable? true, false. If true,  numberOfDesignVariables must be >2
        varyingObjectiveScales % varying objective scales, true, false
        proportionOfNeutralSpace % proprtion of space in 2D set as neutral, 0 <= p < 1
        
        % state initialised based on arguments
        pi1
        pi1Magnitude
        pi2
        pi2Magnitude
        hardConstraintCentres
        hardConstraintRadii
        softConstraintCentres
        softConstraintRadii
        neutralRegionCentres
        neutralRegionRadii
        neutralRegionObjectiveValues = sqrt(8); % maximum distance is between two corners, i.e. sqrt(2^2+2^2)
        discontinuousRegionCentres
        discontinuousRegionRadii
        discontinuousRegionObjectiveValueOffset
        attractorsList
        attractorRegions
        rescaleConstant = 0;
        rescaleMultiplier = 1;
        centreList
        centreRadii
        paretoAngles
        paretoSetIndices % indices in the centreList correcponding to the Pareto set region centres
        rotations
        monte_carlo_samples = 10000;
        pivot_locations
        bracketing_locations_lower
        bracketing_locations_upper
    end
    
    % publically accessible methods of the DBMOPP class
    
    methods
        %--
        
        function obj = DBMOPP(numberOfObjectives, numberOfDesignVariables, numberOfLocalParetoSets, numberOfDominanceResistanceRegions, ...
                numberOfGlobalParetoSets, proportionOfConstrainedSpaceIfChecker, globalParetoSetType, constraintType, ...
                numberOfdiscontinousObjectiveFunctionRegions, variableSolutionDensity, varyingObjectiveScales, proportionOfNeutralSpace,...
                monte_carlo_samples)
            
            % obj = DBMOPP(numberOfObjectives, numberOfDesignVariables, numberOfLocalParetoSets, numberOfDominanceResistanceRegions, ...
            %                numberOfGlobalParetoSets, proportionOfConstrainedSpaceIfChecker, globalParetoSetType, constraintType, ...
            %                numberOfdiscontinuousObjectiveFunctionRegions, variableSolutionDensity, varyingObjectiveScales, proportionOfNeutralSpace)
            %
            % DBMOPP instance constructor -- returns a random instance of a DBMOPP
            % problem, whose characteristics are defined by the arguments
            %
            % setType
            % constraintType -> 0 (no constraint), 1 to 4 hard vertex, centre,
            %       moat and extended checker. 5 to 8 soft vertex, centre, moat
            %       and extended checker.
            % globalParetoSetType -> 0 (duplicate performance), 1 (partially
            %       overlapping performance), 2 (non-intersecting performance)
            
            
            if (numberOfObjectives < 1)
                error('Need numberOfObjectives to be at least 1');
            end
            obj.numberOfObjectives = numberOfObjectives;
            
            if (numberOfDesignVariables < 2)
                error('Need numberOfDesignVariables to be at least 2');
            end
            obj.numberOfDesignVariables = numberOfDesignVariables;
            
            if (numberOfLocalParetoSets < 0)
                error('Need numberOfLocalParetoSets to be at least 1');
            end
            obj.numberOfLocalParetoSets = numberOfLocalParetoSets;
            
            if (numberOfDominanceResistanceRegions < 0)
                error('Need numberOfLocalParetoSets to be at least 1');
            end
            obj.numberOfDominanceResistanceRegions = numberOfDominanceResistanceRegions;
            
            if (numberOfGlobalParetoSets < 1)
                error('Need numberOfGlobalParetoSets to be at least 1');
            end
            obj.numberOfGlobalParetoSets = numberOfGlobalParetoSets;
            
            if (globalParetoSetType ~= 0) && (globalParetoSetType ~= 1) && (globalParetoSetType ~= 2)
                error('globalParetoSetType must be 0, 1 or 2 (duplicate, partially overlapping or non-intersecting performance)');
            end
            obj.globalParetoSetType = globalParetoSetType;
            
            if (constraintType < 0) || (constraintType > 8)
                error('constraintType must be an integer between 0 and 8 inclusive');
            end
            obj.constraintType = constraintType;
            
            if (proportionOfConstrainedSpaceIfChecker < 0) || (proportionOfConstrainedSpaceIfChecker >= 1)
                error('Proportion of the space constrained must be between 0 (inclusive) and 1 (exclusive)');
            end
            if (constraintType ~= 4) && (constraintType ~=8) && (proportionOfConstrainedSpaceIfChecker ~= 0)
                error('Unless constraint type is checker form, proportionOfConstrainedSpaceIfChecker must be 0');
            end
            obj.proportionOfConstrainedSpaceIfChecker = proportionOfConstrainedSpaceIfChecker;
            
            if (numberOfdiscontinousObjectiveFunctionRegions < 0)
                error('numberOfdiscontinousObjectiveFunctionRegions cannot be negative');
            end
            
            if (globalParetoSetType ~= 0) && (numberOfObjectives < 3)
                error('need at least 3 objectives when having partially intersecting or non-overlapping Pareto sets');
            end
            
            obj.numberOfdiscontinousObjectiveFunctionRegions = numberOfdiscontinousObjectiveFunctionRegions;
            
            obj.variableSolutionDensity = variableSolutionDensity;
            obj.varyingObjectiveScales = varyingObjectiveScales;
            
            if (proportionOfNeutralSpace < 0) || (proportionOfNeutralSpace >= 1)
                error('Proportion of the space with a neutral landscape must be between 0 (inclusive) and 1 (exclusive)');
            end
            obj.proportionOfNeutralSpace = proportionOfNeutralSpace;
            
            if ((proportionOfNeutralSpace + proportionOfConstrainedSpaceIfChecker) >=1)
                error('Cannot have a combined proportion of space constrained and neutral greater or equal to 1');
            end
            
            if exist('monte_carlo_samples','var') == false
                obj.monte_carlo_samples = 10000;
            elseif monte_carlo_samples < 1000
                error('Need to set at least 1000 monte carlo samples for approximation');
            else
                obj.monte_carlo_samples = monte_carlo_samples;
            end
            obj.initialise();
        end
        %--
        
        function plotLandscapeForSingleObjective(obj,index,resolution)
            % plotLandscapeForSingleObjective(obj,index,resolution)
            %
            % INPUTS
            %
            % index = index of which objective to plot
            % resolution = mesh size (number of cells on each dimension) to use
            %         when plotting the objective response
            %
            % plots the single objective landscape of the obj instance for the
            % 'index' objective. The optional argument resolution sets the grid
            % on each access (default 500)
            %
            if exist('resolution','var')==false
                resolution = 500;
            end
            
            if resolution < 1
                error('Cannot grid the space with a resolution less than 1');
            end
            if index < 1
                error('Cannot use an index less than 1');
            end
            if index > obj.numberOfObjectives
                error('Cannot use an index greater than the number of objectives in this problem instance');
            end
            
            xy = linspace(-1,1,resolution);
            Z = zeros(resolution,resolution);
            for i=1:resolution
                for j=1:resolution
                    [objective_vector] = obj.evaluate2D([xy(i), xy(j)]);
                    Z(i,j) = objective_vector(index);
                end
            end
            figure;
            surfc(xy,xy,Z');
            shading flat
            view(2)
        end
        %--
        
        function plotParetoSetMembers(obj,resolution)
            % function plotParetoSetMembers(obj,resolution)
            %
            % INPUTS
            %
            % index = index of which objective to plot
            % resolution = mesh size (number of cells on each dimension) to use
            %         when plotting the objective response
            %
            % plots the single objective landscape of the obj instance for the
            % 'index' objective. The optional argument resolution sets the grid
            % on each access (default 500)
            %
            if exist('resolution','var')==false
                resolution = 500;
            end
            
            if resolution < 1
                error('Cannot grid the space with a resolution less than 1');
            end
            
            xy = linspace(-1,1,resolution);
            
            figure;
            hold on;
            for i=1:resolution
                for j=1:resolution
                    if obj.isPareto2D([xy(i), xy(j)])
                        plot(xy(i), xy(j), 'k.');
                    end
                end
            end
            axis([-1 1 -1 1])
            axis square;
        end
        
        %--
        function plotDominanceLandscape(obj,resolution)
            % plotDominanceLandscape(obj,index,resolution)
            %
            % INPUT
            %
            % resolution = mesh size (number of cells on each dimension) to use
            %         when plotting the landscape
            %
            % Plots the instance of this problem
            error('Functionality currently not implemented -- please check the repo often as should be in soon!');
        end
        %--
        function plotProblemInstance(obj)
            % plotProblemInstance(obj)
            %
            % Method plots a visualisation of the problem, showing the location
            % of the attractor points, centres, etc.
            
            figure;
            hold on;
            axis([-1 1 -1 1]);
            axis square
            
            % PLOT ATTRACTOR REGIONS
            
            % plot local Pareto set regions
            for i = 1 : obj.numberOfLocalParetoSets
                %rectangle('Position',[obj.centreList(i,1)-obj.centreRadii(i), obj.centreList(i,2)-obj.centreRadii(i), 2*obj.centreRadii(i), 2*obj.centreRadii(i)],'Curvature',[1, 1],'EdgeColor', [1 0 0])
                C = convhull(obj.attractorRegions{i}.locations);
                fill(obj.attractorRegions{i}.locations(C,1),obj.attractorRegions{i}.locations(C,2),'g');
            end
            
            % plot global Pareto set regions
            for i = obj.numberOfLocalParetoSets+1 : obj.numberOfLocalParetoSets + obj.numberOfGlobalParetoSets
                %rectangle('Position',[obj.centreList(i,1)-obj.centreRadii(i), obj.centreList(i,2)-obj.centreRadii(i), 2*obj.centreRadii(i), 2*obj.centreRadii(i)],'Curvature',[1, 1],'EdgeColor', [0 0 0])
                C = convhull(obj.attractorRegions{i}.locations);
                fill(obj.attractorRegions{i}.locations(C,1),obj.attractorRegions{i}.locations(obj.attractorRegions{i}.convhull,2),'r');
            end
            
            % plot dominance resistance set regions
            for i = obj.numberOfLocalParetoSets + obj.numberOfGlobalParetoSets + 1: obj.numberOfLocalParetoSets + obj.numberOfGlobalParetoSets + obj.numberOfDominanceResistanceRegions
                %rectangle('Position',[obj.centreList(i,1)-obj.centreRadii(i), obj.centreList(i,2)-obj.centreRadii(i), 2*obj.centreRadii(i), 2*obj.centreRadii(i)],'Curvature',[1, 1],'EdgeColor', [0 0 1])
                [n,~] = size(obj.attractorRegions{i}.locations);
                if n>2
                    C = convhull(obj.attractorRegions{i}.locations);
                end
                %obj.attractorRegions{i}.locations
                if n==1
                    plot(obj.attractorRegions{i}.locations(:,1),obj.attractorRegions{i}.locations(:,2),'b.');
                elseif n>2 % enough to draw an area
                    fill(obj.attractorRegions{i}.locations(C,1),obj.attractorRegions{i}.locations(obj.attractorRegions{i}.convhull,2),'b');
                elseif n==2 % just two points, so draw a line
                    plot(obj.attractorRegions{i}.locations(:,1),obj.attractorRegions{i}.locations(:,2),'b-');
                end
            end
            
            % PLOT PENALISED REGIONS
            % plot hard constraint regions
            for i = 1 : length(obj.hardConstraintRadii)
                rectangle('Position',[obj.hardConstraintCentres(i,1)-obj.hardConstraintRadii(i), obj.hardConstraintCentres(i,2)-obj.hardConstraintRadii(i), 2*obj.hardConstraintRadii(i), 2*obj.hardConstraintRadii(i)],'Curvature',[1, 1],'FaceColor', [0 0 0],'EdgeColor', [0 0 0])
            end
            
            % plot soft constraint regions
            for i = 1 : length(obj.softConstraintRadii)
                rectangle('Position',[obj.softConstraintCentres(i,1)-obj.softConstraintRadii(i), obj.softConstraintCentres(i,2)-obj.softConstraintRadii(i), 2*obj.softConstraintRadii(i), 2*obj.softConstraintRadii(i)],'Curvature',[1, 1],'FaceColor', [0.5 0.5 0.5],'EdgeColor', [0.5 0.5 0.5])
            end
            
            % PLOT NEUTRAL REGIONS
            for i = 1 : length(obj.neutralRegionRadii)
                rectangle('Position',[obj.neutralRegionCentres(i,1)-obj.neutralRegionRadii(i), obj.neutralRegionCentres(i,2)-obj.neutralRegionRadii(i), 2*obj.neutralRegionRadii(i), 2*obj.neutralRegionRadii(i)],'Curvature',[1, 1],'FaceColor', [0 1 1],'EdgeColor', [0 1 1])
            end
            
            % PLOT DISCONNECTED PENALTY
            warning("disconnected Pareto penalty regions not yet plotted")
            
            
            % PLOT ATTRACTOR POINTS
            for k = 1:obj.numberOfObjectives
                C = obj.attractorsList{k}.locations;
                plot(C(:,1),C(:,2),'k.');
                text(C(:,1),C(:,2),int2str(k));
            end
            
        end
        
        %--
        
        
        function isPareto = isAParetoSetMember(obj, x, suppressWarning)
            % isPareto = isAParetoSetMember(obj,x,suppressWarning)
            %
            % SHOULD NOT BE USED IN OPTIMISATION PROCESS!!
            %
            % INPUT
            % x = design vector
            % OUTPUT
            % isPareto = true of false depending on if x is an optimal solution
            %
            % Returns true if the design vector passed in is Pareto optimal,
            % otherwise returns false
            if suppressWarning == false
                warning('This function should not be called as part of the optimisation process')
            end
            obj.checkValidLength(x);
            x = obj.get2DVersion(x);
            isPareto = obj.isPareto2D(x);
        end
        
        
        
        function x = getAParetoSetMember(obj, suppressWarning)
            % x = getAParetoSetMember(obj, suppressWarning)
            % SHOULD NOT BE USED IN OPTIMISATION PROCESS!!
            % Returns a random Pareto set member, uniformly from the Pareto set
            if suppressWarning == false
                warning('This function should not be called as part of the optimisation process')
            end
            error('Functionality not currently implemented in this version');
            % while not a legal point obtained
            % get a random Pareto centre
            
            % check for degenerate 2D case
            
            % if centre constraint type used, randomly choose an angle, use
            % corresponding radii from the Pareto set centre list and
            % project and return
            
            % else, generate random point in circle, and check if in convex
            % hull
            
        end
        %--
        
        function [objective_vector, soft_constraint_violation, hard_constraint_violation] = evaluate(obj,x)
            % [objective_vector, soft_constraint_violation, hard_constraint_violation] = evaluate(obj,x)
            %
            % Evalutes the design vector x under this instance of the problem
            %
            % INPUTS
            % x = design vector to evaluate, should adhere to the box
            %   constraints of the problem, i.e. between -1 and +1 on all
            %   dimensions
            %
            % OUTPUTS
            %
            % objective_vector = vector of objective values associated with
            %   x for this DBMOPP problem. If x violates any soft or hard
            %   constraints will be a vector of NaNs
            % soft_constraint_violation = Soft constraint violation. 0 if
            %   no violation, otherwise value increases with distance to the
            %   boundary from the constrain/legal space
            %hard_constraint_violation = Hard constraint violation. Value
            %   of true if there is a violation, false otherwise
            
            obj.checkValidLength(x);
            x = obj.get2DVersion(x);
            [objective_vector, soft_constraint_violation, hard_constraint_violation] = obj.evaluate2D(x);
        end
    end
    
    
    % static helper methods (don't use instance state)
    
    methods(Static)
        % INPUT
        % X1 is a matrix or array, size n x number of design variables
        % x2 is a array, size 1 x number of design variables
        %
        % OUTPUT
        % n-dimensional vector containing the Euclidean distance of each
        % row of X1 to x2
        function d = euclideanDistance(X1,x2)
            d = sqrt(sum((X1-x2).^2,2));
        end
        %--
        % INPUT
        % n = positive integer
        % OUTPUT
        % angles = a vector of n random angles
        function angles = getRandomAngles(n)
            angles = rand(n,1)*2*pi; % arbitrary angles for points on circumference of region
        end
        %--
        % INPUT
        % centres = matrix of regions centres, n by 2
        % radii = vector of radii, 1 by n
        % x = point in Cartesian space, 1 by 2
        % OUTPUT
        % inRegion = true if x is located in any of the circular regions
        %       defined by centres and radii, false otherwise
        % d = array of each of the Euclidean distances from x to each of
        %       the centres
        function [inRegion, d] = inRegion(centres,radii,x)
            inRegion = false;
            d = [];
            if (isempty(centres)==false)
                d = DBMOPP.euclideanDistance(centres,x);
                if sum(d <= radii) > 0
                    inRegion = true;
                end
            end
        end
        %--
        % INPUT
        % x = 2D point to check
        % pivot_location = attractor on boundary of circle
        % location1 = another point on boundary of circle
        % location2 = another point on boundary of circle
        % OUTPUT
        % t = true if x on different side of line defined by pivot_location
        %     and location1, compared to the side of the line defined by
        %     pivot location and location2. If x is also in the circle,
        %     then x is between the two lines if t is true,
        %--
        function t = between_lines_rooted_at_pivot(x, pivot_location, location1, location2)
            t = false;
            d1 = (x(1) - pivot_location(1))*(location1(2) - pivot_location(2)) - (x(2) - pivot_location(2))*(location1(1) - pivot_location(1));
            d2 = (x(1) - pivot_location(1))*(location2(2) - pivot_location(2)) - (x(2) - pivot_location(2))*(location2(1) - pivot_location(1));
            
            if (d1 == 0)
                t = true;
            elseif (d2 == 0)
                t = true;
            elseif sign(d1) ~= sign(d2)
                t = true;
            end
        end
    end
    
    
    % private methods affecting state of object, largely used in the
    % initialisation of a DBMOPP instance
    
    methods(Hidden)
        
        function [objective_vector, soft_constraint_violation, hard_constraint_violation] = evaluate2D(obj,x)
            hard_constraint_violation = obj.getHardConstraintViolation(x);
            if (hard_constraint_violation)
                soft_constraint_violation = 0.0;
                objective_vector = ones(1,obj.numberOfObjectives)*NaN;
                if (obj.constraintType == 3)
                    % need to do an additional check for moat type in case
                    % falls in Pareto set/local front set
                    if (obj.inCOnvexHullOfAtrractorRegion(x))
                        hard_constraint_violation = false;
                        objective_vector = obj.getObjectives(x); % can do this quicker as know region, refactor
                    end
                end
                return;
            end
            
            soft_constraint_violation = obj.getSoftConstraintViolation(x);
            if (soft_constraint_violation > 0)
                hard_constraint_violation = 0.0;
                objective_vector = ones(1,obj.numberOfObjectives)*NaN;
                if (obj.constraintType == 7)
                    % need to do an additional check for moat type in case
                    % falls in Pareto set/local front set
                    if (obj.inCOnvexHullOfAtrractorRegion(x))
                        soft_constraint_violation = 0;
                        objective_vector = obj.getObjectives(x); % can do this quicker as know region, refactor
                    end
                end
                return;
            end
            % neither constraint form is being violated, so calculate the
            % objective vector
            % first check if in neutral region
            inNeutralRegion = DBMOPP.inRegion(obj.neutralRegionCentres, obj.neutralRegionRadii, x);
            % if lies in a neutral region, assign the flat value for
            % objectives
            if inNeutralRegion
                objective_vector = obj.neutralRegionObjectiveValues;
            else
                objective_vector = obj.getObjectives(x);
            end
        end
        
        function isPareto = isPareto2D(obj, x)
            isPareto = obj.is_in_limited_region(x);
            if obj.getHardConstraintViolation(x)
                isPareto = false;
            end
            if obj.getSoftConstraintViolation(x) >0
                isPareto = false;
            end
        end
        
        %--
        % inHull = inConvexHullOfAttractorRegion(obj,x)
        %
        % INPUT
        % x = point to query
        %
        % OUTPUT
        % inHull = true if in convex hull defined by attractor points to induce
        % region of interest (e.g. local/global non-dominated set), false
        % otherwise. Note, does not consider constraints/penalties.
        %
        function inHull = inConvexHullOfAttractorRegion(obj,x)
            obj.checkValidLength(x);
            x = obj.get2DVersion(x); % convert to 2D if problem is in higher dimensions
            
            inHull = false;
            d = DBMOPP.euclideanDistance(obj.centreList,x);
            I = find(d < obj.centreRadii);
            % check if inside any of the circles defining arractor regions
            if isempty(I)==false
                % inside radius of attractor region, need to check if in conv hull
                if inpolygon(x(1),x(2),obj.attractorRegions{I(1)}.locations(obj.attractorRegions{I(1)}.convhull,1), obj.attractorRegions{I(1)}.locations(obj.attractorRegions{I(1)}.convhull,2))
                    inHull = true;
                end
            end
        end
        %--
        function checkValidLength(obj,x)
            if length(x) ~= obj.numberOfDesignVariables
                error(strcat('Number of design variables in the argument does not match that required in the problem instance, was ', int2str(length(x)), ', needs ', int2str(obj.numberOfDesignVariables)));
            end
        end
        %--
        function setUpAttractorCentres(obj)
            n = obj.numberOfLocalParetoSets + obj.numberOfGlobalParetoSets + obj.numberOfDominanceResistanceRegions;
            
            %if regions are grid distributed (most compact), will be sqrt(n) by sqrt(n)
            % need to be 4r apart, and 1r form boundary, and spanning -1 to
            % +1 (i.e. 2), so we can derive max_r from this
            max_r = 1/(2*sqrt(n)+1);
            % now reduce according to the proportion of space to be set
            % aside for constrained/neutral regions
            radius = max_r *(1-(obj.proportionOfNeutralSpace + obj.proportionOfConstrainedSpaceIfChecker));
            
            %place centres of all attractor regions
            %SET AS OPTION ??
            % radius = radius * rand(1,1); % random scaling
            radius = obj.placeRegions(n,radius);
            obj.centreRadii = ones(n,1)*radius;
            %reduce radii if local fronts being used
            if obj.numberOfLocalParetoSets > 0
                globalRadii = radius/2;
                obj.centreRadii(obj.numberOfLocalParetoSets+1:end) = globalRadii;
                w = linspace(1,0.5,obj.numberOfLocalParetoSets+1);
                %linearly decrease local front radii
                obj.centreRadii(1:obj.numberOfLocalParetoSets) = obj.centreRadii(1:obj.numberOfLocalParetoSets).*w(1:obj.numberOfLocalParetoSets)';
            end
            % save indices of Pareto set locations
            obj.paretoSetIndices = obj.numberOfLocalParetoSets+1:obj.numberOfLocalParetoSets+obj.numberOfGlobalParetoSets;
            
        end
        %--
        function [radius] = placeRegions(obj,n,radius)
            effectiveBound = 1-radius;
            threshold = 4*radius;
            obj.centreList = zeros(n,2);
            time_start = tic;
            time_elapsed = 0;
            MAX_ELAPSED = 5; % max seconds before reattempting
            
            % put very first centre in, ensuring set doesn't cross boundary of feasible space
            obj.centreList(1,:) = (rand(1,2)*2*effectiveBound)-(effectiveBound);
            fprintf('Radius: %f\n',radius);
            
            % now assign remaining attractor centre locations
            for i=2:n
                invalid = true;
                while(invalid)
                    c = (rand(1,2)*2*effectiveBound)-(effectiveBound); %random cordinate pair between -(1-radius) and +(1-radius)
                    t = min(DBMOPP.euclideanDistance(obj.centreList(1:i-1,:),c));
                    if ( t > threshold) % centres not equal or closer than 4*radius, so valid to use
                        fprintf('Assigned centre %d\n', i);
                        invalid=false;
                    end
                    time_elapsed = toc(time_start);
                    if (time_elapsed>MAX_ELAPSED)
                        break;
                    end
                end
                obj.centreList(i,:) = c;
            end
            % if not assigned in timeframe, reduce radius and try again
            if (time_elapsed>MAX_ELAPSED)
                fprintf('restarting attractor region placement with smaller radius...\n');
                radius = obj.placeRegions(n, radius*0.95);
            end
        end
        %--
        function placeAttractors(obj)
            
            initial_locations = zeros(obj.numberOfLocalParetoSets+obj.numberOfGlobalParetoSets,2,obj.numberOfObjectives);
            
            % assign attractors per region for local and global fronts
            for i = 1:obj.numberOfLocalParetoSets+obj.numberOfGlobalParetoSets
                locations = repmat(obj.centreList(i,:),obj.numberOfObjectives,1) ...
                    + repmat(obj.centreRadii(i),obj.numberOfObjectives,2).*[cos(obj.paretoAngles + obj.rotations(i)), sin(obj.paretoAngles + obj.rotations(i))];
                obj.attractorRegions{i}.locations = locations;
                obj.attractorRegions{i}.objectiveIndices = 1:obj.numberOfObjectives;
                obj.attractorRegions{i}.centre = obj.centreList(i,:); % some duplication of storage here, remove?
                obj.attractorRegions{i}.radius = obj.centreRadii(i); % some duplication of storage here, remove?
                obj.attractorRegions{i}.convhull = convhull(locations(:,1),locations(:,2));
                for k=1:obj.numberOfObjectives
                    initial_locations(i,:,k) = locations(k,:);
                end
                
            end
            
            % copy across locations used in attractorsList -- we use two
            % structures for effciency later on and for ease of plotting
            % (i.e. obj.attractorsList and obj.attractorRegions)
            for k=1:obj.numberOfObjectives
                obj.attractorsList{k}.locations = zeros(obj.numberOfLocalParetoSets+obj.numberOfGlobalParetoSets,2);
                obj.attractorsList{k}.locations(:,:) = initial_locations(:,:,k);
            end
            
            % now assign dominance resistance regions, which have a subset
            % of attractors per region active
            for i = obj.numberOfLocalParetoSets+obj.numberOfGlobalParetoSets+1 : obj.numberOfLocalParetoSets+obj.numberOfDominanceResistanceRegions+obj.numberOfGlobalParetoSets
                locations = repmat(obj.centreList(i,:),obj.numberOfObjectives,1) ...
                    + repmat(obj.centreRadii(i),obj.numberOfObjectives,2).*[cos(obj.paretoAngles + obj.rotations(i)), sin(obj.paretoAngles + obj.rotations(i))];
                % now flag which objectives to include
                num_to_include = randperm(obj.numberOfObjectives-1);
                num_to_include = num_to_include(1); % number to include
                [~,I] = sort(rand(obj.numberOfObjectives,1)); %random objective indices
                obj.attractorRegions{i}.locations = locations(I(1:num_to_include),:);
                obj.attractorRegions{i}.objectiveIndices = I(1:num_to_include);
                obj.attractorRegions{i}.radius = obj.centreRadii(i);
                for k = 1:num_to_include
                    obj.attractorsList{I(k)}.locations = [obj.attractorsList{I(k)}.locations; locations(I(k),:)];
                end
            end
            
        end
        %--
        % sets up the object instance after parameters are set in the
        % first part of the constructor
        function initialise(obj)
            
            % place attractor centres for regions defining attractor points
            setUpAttractorCentres(obj);
            % set up angles for attractors on regin cicumferences and arbitrary rotations for regions
            obj.paretoAngles = DBMOPP.getRandomAngles(obj.numberOfObjectives); % arbitrary angles for Pareto set
            obj.rotations = rand(length(obj.centreRadii),1)*2*pi;
            % now place attractors
            obj.placeAttractors();
            if obj.globalParetoSetType ~= 0
                obj.placeDisconnectedParetoElements();
            end
            obj.placeDiscontiunitiesNeutralAndCheckerConstraints();
            % set the neural avlue to be the same in all neutral locations
            obj.neutralRegionObjectiveValues = ones(1,obj.numberOfObjectives)*obj.neutralRegionObjectiveValues;
            obj.placeVertexConstraintLocations();
            obj.placeCentreConstraintLocations();
            obj.placeMoatConstraintLocations();
            fprintf('set projection vectors');
            fprintf('set rescaling');
        end
        %--
        function placeDisconnectedParetoElements(obj)
            n = obj.numberOfGlobalParetoSets-1; % number of points to use to set up seperate subregions
            % first get base angles in region of interest on unrotated Pareto set
            
            pivotIndex = randi(obj.numberOfObjectives); % get an attractor at random
            [~, I] = sort(obj.paretoAngles); % sort angles from smallest to largest
            
            if pivotIndex == 1 % if selected as pivot is smallest angle
                offset1Angle =  obj.paretoAngles(I(obj.numberOfObjectives));
            else
                offset1Angle = obj.paretoAngles(I(pivotIndex-1)); % I(pivotIndex)-1;
            end
            if pivotIndex == obj.numberOfObjectives % if selected as pivot is largest
                offset2Angle = obj.paretoAngles(I(1)); %1;
            else
                offset2Angle = obj.paretoAngles(I(pivotIndex+1));
            end
            % angles at offsetIndex1 and offsetIndex2 bracket that at the pivot
            
            %offset1Angle = obj.paretoAngles(offsetIndex1);
            %offset2Angle = obj.paretoAngles(offsetIndex2);
            pivotAngle = obj.paretoAngles(I(pivotIndex));
            
            if (pivotAngle == offset1Angle ) || (pivotAngle == offset2Angle)
                error('angle should not be duplicated');
            end
            
            % now sort, and append each bracketing arractor
            if offset1Angle < offset2Angle % offset angles either side of 2*pi
                % go from 0 to offset1 and offset2 to 2pi
                range_covered = offset1Angle + 2*pi - offset2Angle;
                p1 = offset1Angle /range_covered; % proportion of range in the 0-offset1 angle
                temp = rand(n,1); % drawn n values from Uniform distribution
                r = temp;
                p1 = sum(temp<p1); % get number to draw in range 0-offset1 angle
                r(1:p1,1) = 2*pi + rand(p1,1) * offset1Angle; % for draws less than p1, project to 0-offset1 angle , adding 2*pi as makes sorteding easier later
                r(p1+1:n,1) = rand(n-p1,1) * (2*pi - offset2Angle) + offset2Angle;
                r = sort(r);
                r_angles = zeros(n+1,1);
                r_angles(1,1) = offset2Angle;
                r_angles(n+2,1) = offset1Angle + 2*pi; % adding 2*pi as shifted for sorting
                r_angles(2:n+1,1) = r;
            else
                % go from offset2 to offset1
                r = rand(n,1) * (offset1Angle - offset2Angle) + offset2Angle;
                r = sort(r); %sort from smallest to highest
                r_angles = zeros(n+1,1);
                r_angles(1,1) = offset2Angle;
                r_angles(n+2,1) = offset1Angle;
                r_angles(2:n+1,1) = r;
            end
            % r_angles is now sorted from smallest to largest angle, including
            % offsetIndex2 and offsetIndex1 at either end
            
            
            % now for each Pareto set region, get the corresponding (rotated) locations
            % of the points defining each slice, and save
            
            % preallocate for speed
            obj.pivot_locations = zeros(obj.numberOfLocalParetoSets + obj.numberOfGlobalParetoSets,2);
            obj.bracketing_locations_lower = zeros(obj.numberOfLocalParetoSets + obj.numberOfGlobalParetoSets,2);
            obj.bracketing_locations_upper = zeros(obj.numberOfLocalParetoSets + obj.numberOfGlobalParetoSets,2);
            index = 1;
            
            for i = obj.numberOfLocalParetoSets+1 : obj.numberOfLocalParetoSets + obj.numberOfGlobalParetoSets
                % get the rotated locations for each of the bracketing points
                obj.pivot_locations(i,:) = obj.centreList(i,:) + repmat(obj.centreRadii(i),1,2).*[cos( pivotAngle + obj.rotations(i) ), sin( pivotAngle + obj.rotations(i) )];
                obj.bracketing_locations_lower(i,:) = obj.centreList(i,:) + repmat(obj.centreRadii(i),1,2).*[cos(r_angles(index) + obj.rotations(i)), sin(r_angles(index) + obj.rotations(i))];
                
                if obj.globalParetoSetType == 0
                    error('should not be calling this method with an insatnce with identical Pareto set regions');
                elseif obj.globalParetoSetType == 2 %(non-intersecting performance)
                    obj.bracketing_locations_upper(i,:) = obj.centreList(i,:) + repmat(obj.centreRadii(i),1,2).*[cos(r_angles(index+1) + obj.rotations(i)), sin(r_angles(index+1) + obj.rotations(i))];
                elseif obj.globalParetoSetType == 1 %(partially overlapping performance)
                    if index == obj.numberOfGlobalParetoSets
                        % special case where the relation is inverted
                        obj.bracketing_locations_lower(i,:) = obj.centreList(i,:) + repmat(obj.centreRadii(i),1,2).*[cos(r_angles(2) + obj.rotations(i)), sin(r_angles(2) + obj.rotations(i))];
                        obj.bracketing_locations_upper(i,:) = obj.centreList(i,:) + repmat(obj.centreRadii(i),1,2).*[cos(r_angles(n) + obj.rotations(i)), sin(r_angles(n) + obj.rotations(i))];
                    else
                        obj.bracketing_locations_upper(i,:) = obj.centreList(i,:) + repmat(obj.centreRadii(i),1,2).*[cos(r_angles(index+2) + obj.rotations(i)), sin(r_angles(index+2) + obj.rotations(i))];
                    end
                end
                index = index + 1;
            end
        end
        %--
        function placeVertexConstraintLocations(obj)
            fprintf('Assigning any vertex soft/hard constraint regions\n');
            if (obj.constraintType == 1) || (obj.constraintType == 4)
                % hard or soft vertex
                r = min(obj.centreRadii); % get radius of Pareto Regions -- could store in state to save this calculation?
                % set penalty radius to be no more than half of the raddius
                % of the region associated with an attractor
                penalityRadius = rand(1,1)/2;
                
                [totalToPlace,~] = size(obj.attractorsList{1});
                for i=2:length(obj.attractorsList)
                    [t,~] = size(obj.attractorsList{i});
                    totalToPlace = totalToPlace + t;
                end
                % preallocate for speed
                centres = zeros(totalToPlace,2);
                radii = zeros(totalToPlace,1);
                k=1;
                for i=1:length(obj.attractorRegions)
                    for j=1:length(obj.attractorRegions{i}.objectiveIndices)
                        centres(k,:) = obj.attractorRegions{i}.locations(j,:);
                        radii(k) = obj.attractorRegions{i}.radius * penalityRadius;
                        k = k+1;
                    end
                end
                % now assign to hard or soft
                if (obj.constraintType == 1)
                    obj.hardConstraintCentres = centres;
                    obj.hardConstraintRadii = radii;
                else
                    obj.softConstraintCentres = centres;
                    obj.softConstraintRadii = radii;
                end
                
            end
        end
        %--
        function placeCentreConstraintLocations(obj)
            fprintf('Assigning any centre soft/hard constraint regions\n');
            if (obj.constraintType == 2)
                obj.hardConstraintCentres = obj.centreList;
                obj.hardConstraintRadii = obj.centreRadii;
            elseif (obj.constraintType == 5)
                obj.softConstraintCentres = obj.centreList;
                obj.softConstraintRadii = obj.centreRadii;
            end
        end
        %--
        function placeMoatConstraintLocations(obj)
            fprintf('Assigning any moat soft/hard constraint regions\n');
            if (obj.constraintType == 3) || (obj.constraintType == 5)
                % need to consider how wide to set the moat, going to set as
                % random between 0 and radius of the max region size (will
                % ensure no overlap with any Pareto regions)
                
                r = rand(1,1)+1;
                if (obj.constraintType == 3)
                    obj.hardConstraintCentres = obj.centreList;
                    obj.hardConstraintRadii = obj.centreRadii*r;
                elseif (obj.constraintType == 6)
                    obj.softConstraintCentres = obj.centreList;
                    obj.softConstraintRadii = obj.centreRadii*r;
                end
                
            end
        end
        %--
        function placeDiscontiunitiesNeutralAndCheckerConstraints(obj)
            fprintf('Assigning any checker soft/hard constraint regions and neutral regions\n');
            if ((obj.proportionOfConstrainedSpaceIfChecker + obj.proportionOfNeutralSpace) > 0)
                
                % at this point the attractor regions are set, but we need to
                % place the (non-attractor-intersecting) regions for netrality and
                % extended check-type constrains, and also disconuity in the
                % functions
                
                % proportionOfConstrainedSpaceIfChecker proportionOfNeutralSpace
                
                % Monte Carlo samples of space to limit neutral and constrained space amount
                
                
                S = (rand(obj.monte_carlo_samples,2)*2) -1; %samples in domain to track total regions
                
                % first remove any samples falling in attractor regions for
                % computational efficiency
                for i=1:length(obj.centreRadii)
                    d = DBMOPP.euclideanDistance(S,obj.centreList(i,:));
                    S((d<=obj.centreRadii(i)),:)=[]; % remove samples falling within region radius of
                end
                if (length(S) < obj.monte_carlo_samples*(obj.proportionOfConstrainedSpaceIfChecker + obj.proportionOfNeutralSpace) )
                    error('Not enough space outside of attractor regions to match requirement of constrained+neural space');
                end
                % Now iteratively place a centre, check legality, and update
                % proportion of space covered as estimated using the MC samples
                % falling inside the neutral/penality region
                % Note, by definition, all samples in S are legal centres
                % outside of attractor regions, and are randomly ordered
                % so whill just select from these to speed up the process
                
                % declare outside, as will be used by neutral space
                % allocation code
                constrainedCentres = [];
                constrainedRadii = [];
                % do constrained space first
                if (obj.proportionOfConstrainedSpaceIfChecker > 0)
                    [constrainedCentres,constrainedRadii,S] = obj.setNotAttractorRegionsAsProportionOfSpace(S,obj.proportionOfConstrainedSpaceIfChecker,[],[]);
                    
                    if obj.constraintType == 4 % hard checker
                        obj.hardConstraintCentres = constrainedCentres;
                        obj.hardConstraintRadii = constrainedRadii;
                    else  % soft checker
                        if (obj.constraintType ~= 8)
                            error('constraintType shuld be 8 to reach here');
                        end
                        obj.softConstraintCentres = constrainedCentres;
                        obj.softConstraintRadii = constrainedRadii;
                    end
                end
                
                %now do neutral space
                if (obj.proportionOfNeutralSpace > 0 )
                    [centres,radii] = obj.setNotAttractorRegionsAsProportionOfSpace(S,obj.proportionOfNeutralSpace,constrainedCentres,constrainedRadii);
                    obj.neutralRegionCentres = centres;
                    obj.neutralRegionRadii = radii;
                end
                % disconunity placement
                fprintf('Check discontinuity');
            end
        end
        %--
        function [constrainedCentres,constrainedRadii,S] = setNotAttractorRegionsAsProportionOfSpace(obj,S,proportionToAttain,otherCentres,otherRadii)
            allocated = 0;
            constrainedCentres = [];
            constrainedRadii = [];
            % iteratively add a constrained area
            while allocated < proportionToAttain
                constrainedCentres = [constrainedCentres; S(end,:)];
                % distance to the centres of the regions, and those
                % allocated in the 'otherCentres' already allocated constrained areas
                d = DBMOPP.euclideanDistance([obj.centreList; otherCentres], S(end,:));
                d = d - [obj.centreRadii; otherRadii];
                d = min(d);
                % SANITY CHECK
                if (d<=0)
                    error('This state should not occur');
                end
                
                % d now holds the maximum radius that a penality
                % region centred on S(end,:) may have, before it
                % conflicts with an arractor region
                % Now calculate the radius to use to cover remain
                % proportion desired
                %
                % pi * max_r ^2 = proportionStillToCover
                % c_r = sqrt(proportionStillToCover/pi)
                
                c_r = sqrt((proportionToAttain-allocated)/pi);
                % set max radius to use to be the minimum of the
                % radius infringing the nearest region, and that
                % required to reach the target covered
                r = rand(1,1)*min(d,c_r);
                
                constrainedRadii = [constrainedRadii; r];
                S(end,:) = []; % remove allocated point from S
                % now remove any other covered points
                d = DBMOPP.euclideanDistance(S,constrainedCentres(end,:));
                I = d<=r;
                S(I,:)=[];
                coveredCount = 1+sum(I); % count the number removed from S
                
                % update the approximation of the proportion of
                % space covered with the fraction of space
                % *exclusively* covered with the new region
                allocated = allocated + coveredCount/obj.monte_carlo_samples;
            end
            % now constrainedCentres and constrainedRadii cover
            % collectively proportionOfConstrainedSpaceIfChecker of
            % the design space (in 2D) and do not cnflict with the
            % attractor regions.
            
        end
        %--
        function h = getHardConstraintViolation(obj,x)
            inHardConstraintRegion = DBMOPP.inRegion(obj.hardConstraintCentres, obj.hardConstraintRadii, x);
            if inHardConstraintRegion
                h=1;
            else
                h=0;
            end
        end
        %--
        function s = getSoftConstraintViolation(obj,x)
            [inSoftConstraintRegion, d] = DBMOPP.inRegion(obj.softConstraintCentres, obj.softConstraintRadii, x);
            
            if inSoftConstraintRegion
                k = sum(d <= obj.softConstraintRadii);
                if (k > 0)
                    c = d - obj.softConstraintRadii; % get a measure of how far away from the boundary of the region a point is
                    c = c.*k; % only count those *inside* the constraint region
                    s = max(c);
                else
                    s = 0.0;
                end
            else
                s = 0.0;
            end
        end
        %--
        function x = get2DVersion(obj,x)
            if length(x)>2
                x = [(x*obj.pi1)/obj.pi1Magnitude, (x*obj.pi2)/obj.pi2Magnitude];
            end
        end
        %--
        function y = getMinimumDistancesToAttractors(obj,x)
            y = zeros(1, obj.numberOfObjectives);
            for i = 1:obj.numberOfObjectives
                d = DBMOPP.euclideanDistance(obj.attractorsList{i}.locations,x);
                y(i) = min(d);
            end
            % apply and objectives rescaling
            y = y .* obj.rescaleMultiplier;
            y = y + obj.rescaleConstant;
        end
        %--
        function y = getObjectives(obj,x)
            if (obj.globalParetoSetType == 0)
                % simple case where just need distances from attractors
                y = obj.getMinimumDistancesToAttractors(x);
            else
                y = obj.getMinimumDistancesToAttractorsOverlapOrDiscontinuousForm(x);
            end
            y = obj.updateWithDiscontinuity(x,y);
            y = obj.updateWithNeutrality(x,y);
        end
        %--
        function y = getMinimumDistancesToAttractorsOverlapOrDiscontinuousForm(obj,x)
            y = obj.getMinimumDistancesToAttractors(x);
            [inParetoRegion, inHull, index] = obj.is_in_limited_region(x);
            if inHull
                if ~inParetoRegion
                    y = y + obj.centreRadii(index); % penalise by radius
                end
            end
        end
        %--
        function [inParetoRegion, inHull, index] = is_in_limited_region(obj,x)
            inParetoRegion = false;
            inHull = false;
            d = DBMOPP.euclideanDistance(obj.centreList,x);
            I = find(d < obj.centreRadii);
            % check if inside any of the circles defining arractor regions
            index = -1;
            if isempty(I)==false
                % check if in Global Pareto region rather than e.g. local set
                if (I(1) > obj.numberOfLocalParetoSets) && (I(1) <= obj.numberOfLocalParetoSets + obj.numberOfGlobalParetoSets)
                    % 'inside radius of global attractor region, need to check if in conv hull'
                    if inpolygon(x(1),x(2),obj.attractorRegions{I(1)}.locations(obj.attractorRegions{I(1)}.convhull,1), obj.attractorRegions{I(1)}.locations(obj.attractorRegions{I(1)}.convhull,2))
                        inHull = true;
                    end
                end
            end
            
            if obj.globalParetoSetType == 0 % identical performance
                inParetoRegion = inHull;
            else
                if inHull
                    index = I(1);
                    inParetoRegion = DBMOPP.between_lines_rooted_at_pivot(x, obj.pivot_locations(I(1),:), obj.bracketing_locations_lower(I(1),:), obj.bracketing_locations_upper(I(1),:) );
                    if obj.globalParetoSetType == 1 % partially overlapping performance
                        if I(1) == obj.numberOfLocalParetoSets + obj.numberOfGlobalParetoSets
                            inParetoRegion = ~inParetoRegion; % special case where last region is split at the two sides
                        end
                    end
                end
            end
            
        end
        %--
        function y = updateWithDiscontinuity(obj,x,y)
            if (isempty(obj.discontinuousRegionCentres)==false)
                d = DBMOPP.euclideanDistance(obj.discontinuousRegionCentres,x);
                v = d < obj.discontinuousRegionRadii;
                d = d.*v;
                if sum(d)>0
                    [~,i] = min(d); % choose centre closest to
                    y = y + obj.discontinuousRegionObjectiveValueOffset(i,:);
                end
            end
        end
        %--
        function y = updateWithNeutrality(obj,x,y)
            if (isempty(obj.neutralRegionCentres)==false)
                d = DBMOPP.euclideanDistance(obj.neutralRegionCentres,x);
                v = d < obj.neutralRegionRadii;
                d = d.*v;
                if sum(d)>0
                    [~,i] = min(d); %choose centre closest to
                    y = obj.neutralRegionObjectiveValues(i,:);
                end
            end
        end
    end
    
end

