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
            %                numberOfdiscontinuousObjectiveFunctionRegions, variableSolutionDensity, varyingObjectiveScales, proportionOfNeutralSpace,...
            %                monte_carlo_samples)
            %
            % DBMOPP instance constructor -- returns a random instance of a DBMOPP
            % problem, whose characteristics are defined by the arguments
            %
            % INPUTS
            %
            % numberOfObjectives = number of objectives in instance
            % numberOfDesignVariables = number of design variables (minimum 2)
            % numberOfLocalParetoSets = number of local Pareto sets, minimum 0
            % numberOfDominanceResistanceRegions = number of dominance
            %     resistance regions, minimum 0
            % numberOfGlobalParetoSets = number of Pareto sets (disconnected),
            %     minimum 1
            % proportionOfConstrainedSpaceIfChecker = if extended checker constraint
            %     type is used, proportion of 2D design space to be of this
            %     type
            % globalParetoSetType = 0 (duplicate performance), 1 (partially
            %       overlapping performance), 2 (non-intersecting performance)
            % constraintType = 0 (no constraint), 1 to 4 hard vertex, centre,
            %       moat and extended checker. 5 to 8 soft vertex, centre, moat
            %       and extended checker
            % numberOfdiscontinuousObjectiveFunctionRegions = number of
            %     regions to apply whose boundaries cause discontinuities
            %     in objective functions
            % variableSolutionDensity = Should solution density vary in
            %     mapping down to each of the two visualised dimensions?
            %     (true or false)
            % varyingObjectiveScales = are objective scale varied? (true or
            %     false)
            % proportionOfNeutralSpace = proportion of 2D space to make
            %     neutral
            % monte_carlo_samples = default 10000
            
            
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
            surfc(xy,xy,Z'); %the x-coordinates of the vertices corresponding to column indices of Z and the y-coordinates corresponding to row indices of Z, so transposing
            view(2)
            shading flat
            axis square
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
        
        function [neutral_areas, dominated, Y, destination, dominating_neighbours, offset, basins] = plotDominanceLandscape(obj,resolution,moore_neighbourhood)
            % [neutral_areas, dominated, Y, destination, dominating_neighbours, offset] = plotDominanceLandscape(obj,resolution,moore_neighbourhood)
            %
            % INPUTS
            %
            % index = index of which objective to plot
            % resolution = mesh size (number of cells on each dimension) to use
            %         when plotting the objective response
            % moore_neighbourhood = type of neighbourhood used, if true then
            %         Moore neighbourhood (default), if false Von Neumann
            %         nieghbourhood
            %
            % OUTPUTS
            %
            % neural_areas = resolution by resolution matrix. contigious
            %       dominance neutral araes have same positive integer value. 
            %       dominated located have a value of -1
            % dominated = Boolean resolution by resolution matrix. true is
            %       corresponding location is dominated.
            % Y = number of objectives by resolution by resolution matrix
            %       holding objective vector for each mesh location
            % destination = resolution by resolution cell matrix, holding
            %       the list of distinct neutral areas reached by all
            %       downhill dominance walks commences at the cell (see the
            %       neutral_areas matrix for mapping)
            % offset = negihbourhood mapping matrix 
            % basins = matrix of basin memberships (resolution by 
            %       resolution). A value of 0 at basins(i,j) denotes the
            %       corresponding location is dominance neutral. A value of
            %       1 denotes that dominance paths rooted at the cell lead
            %       to more than one distinct neutral region. A value
            %       between 0.25 and 0.75 denotes that all dominance paths
            %       starting from this location end in the same dominance
            %       neutral regions -- all members of the same basin having
            %       the same value in this range.
            %
            % plots the single objective landscape of the obj instance for the
            % 'index' objective. The optional argument resolution sets the grid
            % on each access (default 500)
            %
            if exist('resolution','var')==false
                resolution = 500;
            end
            if exist('moore_neighbourhood','var')==false
                moore_neighbourhood = true;
            end
            if resolution < 1
                error('Cannot grid the space with a resolution less than 1');
            end
            
            xy = linspace(-1,1,resolution);
            Y = zeros(obj.numberOfObjectives,resolution,resolution);
            % evaluate alll locations
            for i=1:resolution
                for j=1:resolution
                    [objective_vector] = obj.evaluate2D([xy(i), xy(j)]);
                    Y(:,i,j) = objective_vector;
                end
            end
            [neutral_areas, dominated, destination, dominating_neighbours, offset, basins] = DBMOPP.plotDominanceLandscapeFromMatrix(Y,xy,xy,moore_neighbourhood);
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
                C = convhull(obj.attractorRegions{i}.locations);
                fill(obj.attractorRegions{i}.locations(C,1),obj.attractorRegions{i}.locations(C,2),'g');
            end
            
            % plot global Pareto set regions
            for i = obj.numberOfLocalParetoSets+1 : obj.numberOfLocalParetoSets + obj.numberOfGlobalParetoSets
                [n,~] = size(obj.attractorRegions{i}.locations);
                if n>2
                    C = convhull(obj.attractorRegions{i}.locations);
                    fill(obj.attractorRegions{i}.locations(C,1),obj.attractorRegions{i}.locations(obj.attractorRegions{i}.convhull,2),'r');
                else % just two points, so draw a line
                    plot(obj.attractorRegions{i}.locations(:,1),obj.attractorRegions{i}.locations(:,2),'r-');
                end
            end
            
            % plot dominance resistance set regions
            for i = obj.numberOfLocalParetoSets + obj.numberOfGlobalParetoSets + 1: obj.numberOfLocalParetoSets + obj.numberOfGlobalParetoSets + obj.numberOfDominanceResistanceRegions
                [n,~] = size(obj.attractorRegions{i}.locations);
                if n>2
                    C = convhull(obj.attractorRegions{i}.locations);
                end
                if n==1
                    % just draw the point
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
        
        function [x, corresponding2Dpoint] = getAParetoSetMember(obj, suppressWarning)
            % [x, corresponding2Dpoint] = getAParetoSetMember(obj, suppressWarning)
            % SHOULD NOT BE USED IN OPTIMISATION PROCESS!!
            % Returns a random Pareto set member, uniformly from the Pareto
            % set, and the point in 2D it maps to
            
            if suppressWarning == false
                warning('This function should not be called as part of the optimisation process')
            end
            % while not a legal point obtained
            % get a random Pareto centre
            
            invalid = true;
            x = [];
            while invalid
                k = randi(obj.numberOfGlobalParetoSets) + obj.numberOfLocalParetoSets;
                
                % check for degenerate 2D case
                if (obj.constraintType == 2) || (obj.constraintType == 6)
                    % if centre constraint type used, randomly choose an angle, use
                    % corresponding radii from the Pareto set centre list and
                    % project and return
                    angle = rand() * 2.0 * pi;
                    x = obj.centreList(k,:) + [obj.centreRadii(k) * cos(angle), obj.centreRadii(k)*sin(angle)];
                    invalid = false;
                else
                    % else, generate random point in circle,
                    r = obj.centreRadii(k) * sqrt(rand());
                    angle = rand() * 2.0 * pi;
                    x = obj.centreList(k,:) + [r * cos(angle), r*sin(angle)];
                    if (obj.isPareto2D(x)) % check if sample is Pareto optimal
                        invalid = false;
                    end
                end
            end
            % project to higher number of dimensions if required
            if obj.numberOfDesignVariables > 2
                % design space is large than 2D, so need to randomly select
                % a location in this higher dimensional space which maps to
                % this Pareto optimal location
                corresponding2Dpoint = x;
                x = DBMOPP.get_vectors_mapping_to_location(x,obj.pi1,obj.pi1Magnitude,obj.pi2,obj.pi2Magnitude,obj.numberOfDesignVariables);
            else
                corresponding2Dpoint = x;
            end
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
        function [neutral_areas, dominated, destination, dominating_neighbours, offset, basins] = plotDominanceLandscapeFromMatrix(Y,x,y,moore_neighbourhood)
            % [neutral_areas, dominated, destination, dominating_neighbours, offset] = plotDominanceLandscapeFromMatrix(Y,x,y,moore_neighbourhood)
            %
            % INPUTS
            %
            % Y = number of objectives by resolution by resolution matrix
            %       holding objective vector for each mesh location
            % x = resolution by 1 array of ordered x locations of samples
            %       (e.g. x = linspace(-1, 1, resolution) )
            % y = resolution by 1 array of ordered x locations of samples
            %       (e.g. y = linspace(-1, 1, resolution) )
            % moore_neighbourhood = type of neighbourhood used, if true then
            %         Moore neighbourhood, if false Von Neumann
            %         nieghbourhood used
            %
            % OUTPUTS
            %
            % neural_areas = resolution by resolution matrix. contigious
            %       dominance neutral araes have same positive integer value. 
            %       dominated located have a value of -1
            % dominated = Boolean resolution by resolution matrix. true is
            %       corresponding location is dominated.
            % destination = resolution by resolution cell matrix, holding
            %       the list of distinct neutral areas reached by all
            %       downhill dominance walks commences at the cell (see the
            %       neutral_areas matrix for mapping)
            % offset = neighbourhood mapping matrix 
            % basins = matrix of basin memberships (resolution by 
            %       resolution). A value of 0 at basins(i,j) denotes the
            %       corresponding location is dominance neutral. A value of
            %       1 denotes that dominance paths rooted at the cell lead
            %       to more than one distinct neutral region. A value
            %       between 0.25 and 0.75 denotes that all dominance paths
            %       starting from this location end in the same dominance
            %       neutral regions -- all members of the same basin having
            %       the same value in this range.
            %
            [basins, neutral_areas, dominated, destination, dominating_neighbours, offset] = DBMOPP.getDominanceLandscapeBasinsFromMatrix(Y,x,y,moore_neighbourhood);
            
            % in basins all neutral areas have value 1
            % all locations leading to multiple local optima have value 1
            % all locations leading to the same single optima have the same
            % value
            figure;
            imagesc(x,y,basins');
            colormap(gray)
            set(gca,'YDir','normal')
            view(2)
            axis square
        end
        
        function [basins, neutral_areas, dominated, destination, dominating_neighbours, offset] = getDominanceLandscapeBasinsFromMatrix(Y,x,y,moore_neighbourhood)
            % [basins, neutral_areas, dominated, destination, dominating_neighbours, offset] = getDominanceLandscapeBasinsFromMatrix(Y,x,y,moore_neighbourhood)
            %
            % INPUTS
            %
            % Y = number of objectives by resolution by resolution matrix
            %       holding objective vector for each mesh location
            % x = resolution by 1 array of ordered x locations of samples
            %       (e.g. x = linspace(-1, 1, resolution) )
            % y = resolution by 1 array of ordered x locations of samples
            %       (e.g. y = linspace(-1, 1, resolution) )
            % moore_neighbourhood = type of neighbourhood used, if true then
            %         Moore neighbourhood, if false Von Neumann
            %         nieghbourhood used
            %
            % OUTPUTS
            %
            % basins = matrix of basin memberships (resolution by 
            %       resolution). A value of 0 at basins(i,j) denotes the
            %       corresponding location is dominance neutral. A value of
            %       1 denotes that dominance paths rooted at the cell lead
            %       to more than one distinct neutral region. A value
            %       between 0.25 and 0.75 denotes that all dominance paths
            %       starting from this location end in the same dominance
            %       neutral regions -- all members of the same basin having
            %       the same value in this range.
            % neural_areas = resolution by resolution matrix. contigious
            %       dominance neutral araes have same positive integer value. 
            %       dominated located have a value of -1
            % dominated = Boolean resolution by resolution matrix. true is
            %       corresponding location is dominated.
            % destination = resolution by resolution cell matrix, holding
            %       the list of distinct neutral areas reached by all
            %       downhill dominance walks commences at the cell (see the
            %       neutral_areas matrix for mapping)
            % offset = neighbourhood mapping matrix 
            %
            % plots dominance landscae given arguments
            [numObjectives,resolution,r] = size(Y);  
            if (resolution ~= r)
                error('Second and third dimension of Y must be the same size');
            end
            if (numObjectives < 2)
                error('Must have at least two objectives')
            end
            if length(x) ~= resolution
               error('must be as many x grid labels as elements'); 
            end
            if length(y) ~= resolution
               error('must be as many y grid labels as elements'); 
            end
            
            if moore_neighbourhood == true
                dominating_neighbours = false(resolution,resolution,8);
                neutral_neighbours = false(resolution,resolution,8);
            else
                dominating_neighbours = false(resolution,resolution,4);
                neutral_neighbours = false(resolution,resolution,4);
            end
            dominated = true(resolution,resolution);
            
            % array holding neighbourhood directions
            offset = [1 0; -1 0; 0 1; 0 -1; 1 1; 1 -1; -1 1; -1 -1];
            % determine those which are not dominated
            for i=1:resolution
                for j=1:resolution
                    [dominating_neighbours(i,j,:), neutral_neighbours(i,j,:), dominated(i,j)] = DBMOPP.identify_dominating_neighbours(Y,i,j,resolution,moore_neighbourhood,offset);
                end
            end
            % dominating_neighbours now holds location of neighbours which
            % dominate, and dominated holds whether a particular location
            % is dominated by any nieghbour.
            neutral_areas = ones(resolution,resolution) * -1;
            neutral_index = 1;
            % label contiguous neutral areas
            for i=1:resolution
                for j=1:resolution
                    % if not-dominated, and not yet assigned to a neutral
                    % area
                    if (dominated(i,j) == false) && (neutral_areas(i,j) < 0)
                        neutral_areas = DBMOPP.identify_neutral_area_members(i,j,neutral_areas,neutral_index,dominated,neutral_neighbours,moore_neighbourhood,offset);
                        neutral_index = neutral_index + 1;
                    end
                end
            end
            % identify basins of attraction
            basins = ones(resolution,resolution)*-1;
            basins(neutral_areas>0) = 0;
            processed = false(resolution,resolution);
            number_distinct_neutral_regions = max(max(neutral_areas));
            destination = cell(resolution,resolution);
            for i=1:resolution
                for j=1:resolution
                    if ~processed(i,j) % if not yet processed
                        [processed,destination,~] = DBMOPP.update_destinations(i,j,processed,destination,dominating_neighbours,neutral_areas,offset);
                    end
                end
            end
            % now fill value in basins
            for i=1:resolution
                for j=1:resolution
                    if (neutral_areas(i,j) > 0)
                        basins(i,j) = 0;
                    else
                        if length(destination{i,j})==1
                            % put between 0.25 and 0.75, graded by which basin it leads to
                            basins(i,j) = 0.25 + destination{i,j}/( 2*number_distinct_neutral_regions );
                        else
                            basins(i,j) = 1;
                        end
                    end
                end
            end
        end
        
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
        % INPUT
        % centres = matrix of regions centres, n by 2
        % radii = vector of radii, 1 by n
        % x = point in Cartesian space, 1 by 2
        % OUTPUT
        % inRegion = true if x is located in any of the circular regions
        %       defined by centres and radii, false otherwise
        % d = array of each of the Euclidean distances from x to each of
        %       the centres
        function [inRegion, d] = inRegionExcludingBoundary(centres,radii,x)
            inRegion = false;
            d = [];
            if (isempty(centres)==false)
                d = DBMOPP.euclideanDistance(centres,x);
                if sum(d < radii) > 0
                    inRegion = true;
                end
            end
        end
        
        function  plotDominatedWalks(x_index,y_index,dominating_neighbours, offset, neutral_areas)
            [n1,n2,~] = size(dominating_neighbours);
            if (n1~=n2)
                error('input does not match expected, should be a square grid');
            end
            if (x_index > n1) || (y_index > n1)
                error('cannot use an index larger than the grid size');
            end
            paths = ones(n1,n2);
            paths(x_index,y_index) = 0; % highlight starting point
            paths = DBMOPP.recursive_fill_paths(paths,x_index,y_index,dominating_neighbours, offset, false(n1,n2), neutral_areas);
            
            xy = linspace(-1,1,n1);
            
            figure;
            imagesc(xy,xy,paths');
            colormap(gray)
            set(gca,'YDir','normal')
            view(2)
            axis square
        end
    end
    
    methods(Static,Hidden)
        function [processed,destination,destination_list] = update_destinations(i,j,processed,destination,dominating_neighbours,neutral_areas,offset)
            destination_list = [];
            if sum(dominating_neighbours(i,j,:)) > 0 % the cell at i and j has dominating neighbours
                for k = 1 : length(dominating_neighbours(i,j,:))
                    if dominating_neighbours(i,j,k) % query each dominating neighbour
                        if processed(i + offset(k,1),j + offset(k,2)) % if already processed
                            % already know where neighbour leads, so
                            % update destination list with neutral areas
                            % led to by all dominating walks from this
                            % neighbouring dominating location
                            individual_destination_list = destination{i + offset(k,1), j + offset(k,2)};
                        else
                            % neighbour k is dominating, so take all paths downhill
                            [processed,destination,individual_destination_list] = DBMOPP.update_destinations(i + offset(k,1),j + offset(k,2),processed,destination,dominating_neighbours,neutral_areas,offset);
                        end
                        destination_list = unique([destination_list individual_destination_list]);
                    end
                end
                processed(i,j) = true;
                destination{i,j} = destination_list; % record all neutral areas reachable from this point
            else % no dominating neighbours, so at a dominance neutral point
                processed(i,j) = true;
                destination_list = neutral_areas(i,j); % use value at this location, as reached a neutral point
                destination{i,j} = destination_list;
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
        
        function z = get_vectors_mapping_to_location(x,pi1,pi1_mag,pi2,pi2_mag,dim)
            
            z = zeros(1,dim);
            z = DBMOPP.process_dimensions(z,x(1),pi1,pi1_mag);
            z = DBMOPP.process_dimensions(z,x(2),pi2,pi2_mag);
        end
        
        function z = process_dimensions(z,x,pi,pi_mag)
            if pi_mag == 1
                % special case
                z(pi) = x;
            else
                % map value from [-1 +1] to [0 1]
                x = ((x+1)/2)*pi_mag;
                s = DBMOPP.unit_hypercube_simplex_sample(pi_mag,x);
                % map s back to [-1 +1]
                s = (s*2)-1;
                z(pi) = s;
            end
        end
        
        function [paths, path_taken] = recursive_fill_paths(paths, i, j, dominating_neighbours, offset, path_taken, neutral_areas)
            if sum(dominating_neighbours(i,j,:)) > 0 % the cell at i and j has dominating neighbours
                for k = 1 : length(dominating_neighbours(i,j,:))
                    if dominating_neighbours(i,j,k) % query each dominating neighbour
                        if ~path_taken(i+offset(k,1),j+offset(k,2)) % if
                            % at a neighbour which dominates the current
                            % point, which hasn't yet been passed through,
                            % so colour the cell, set as seen and process
                            % all of its dominating neighbour paths
                            paths(i+offset(k,1),j+offset(k,2)) = 0.5;
                            path_taken(i+offset(k,1),j+offset(k,2)) = true;
                            [paths, path_taken] = DBMOPP.recursive_fill_paths(paths,i+offset(k,1),j+offset(k,2),dominating_neighbours, offset,path_taken, neutral_areas);
                        end
                    end
                end
            end
        end
        
        function [dominating_neighbours, neutral_neighbours, dominated] = identify_dominating_neighbours(Y,i,j,resolution,moore_neighbourhood,offset)
            if moore_neighbourhood
                n = 8;
            else
                n = 4;
            end
            dominating_neighbours = false(n,1);
            neutral_neighbours = false(n,1);
            
            for k=1:n
                if (i + offset(k,1) > 0) && (i + offset(k,1) <= resolution) && (j + offset(k,2) > 0) && (j + offset(k,2) <= resolution)
                    [dominating_neighbours(k),neutral_neighbours(k)] = DBMOPP.vector_is_dominated_or_neutral(Y(:,i,j),Y(:,i + offset(k,1),j + offset(k,2)));
                end
            end
            dominated = any(dominating_neighbours);
        end
        
        function neutral_areas = identify_neutral_area_members(i,j,neutral_areas,neutral_index,dominated,neutral_neighbours,moore_neighbourhood,offset)
            % the location at i, j is not dominated, so want to identify
            % all members of contiguous non-dominated area and label
            % accordingly.
            neutral_areas(i,j) = neutral_index;
            for k = 1:length(neutral_neighbours(i,j,:))
                if (neutral_neighbours(i,j,k)) % if neighbour is neutral with respect to current location
                    if neutral_areas(i + offset(k,1), j + offset(k,2)) < 0 % if neighbour not already labelled
                        if (dominated(i + offset(k,1), j + offset(k,2)) == false) % and is not dominated
                            % recursively go through
                            neutral_areas = DBMOPP.identify_neutral_area_members(i + offset(k,1), j + offset(k,2),neutral_areas,neutral_index,dominated,neutral_neighbours,moore_neighbourhood,offset);
                        end
                    end
                end
            end
            
        end
        
        function d = vector_dominates(y1,y2)
            d = (sum(y1 <= y2) == length(y1)) && (sum(y1 < y2) > 0);
        end
        
        function [d, n] = vector_is_dominated_or_neutral(y1,y2)
            % d is true if y2 dominates y1
            % n is true if y1 and y2 are incomparable under the dominates
            % relation
            d = DBMOPP.vector_dominates(y2,y1);
            if d == false
                n = ~DBMOPP.vector_dominates(y1,y2);
            else
                n = false;
            end
        end
        
        function X = unit_hypercube_simplex_sample(dim,sum_value,number_of_points)
            
            % function X = unit_hypercube_simplex_sample(number_of_points, dim, sum_value)
            %
            % INPUTS
            %
            % dim = dimensions
            % sum_value = value that each vector should sum to
            % number_of_points = number of dim-dimensional points to sample (output in
            %       X)
            %
            % OUTPUT
            %
            % X = number_of_points by dim matrix of uniform samples in unit cube from
            %       simplex which sums to sum_value
            %
            % Returns X which contains uniform samples from the simplex summing to
            % sum_value, which lies in the unit hypercube
            
            
            if dim < sum_value
                error('It is impossible to attain sum_value, given variables are in unit range. The sum value given is greater than the number of dimensions');
            end
            if sum_value < 0
                error('Cannot have negative sum value, as variables are in unit range ')
            end
            if exist('number_of_points','var')==false
                number_of_points = 1;
            end
            if number_of_points < 1
                error('must require a positive number of points')
            end
            
            X = exprnd(ones(number_of_points,dim));
            S = sum(X,2);
            if sum_value == 1
                X = X./repmat(S,1,dim);
            elseif sum_value < 1
                X = (X./repmat(S,1,dim))*sum_value;
            else
                if sum_value < dim/2
                    % rejection sampling
                    X = DBMOPP.recalibrate(X,number_of_points,S,sum_value,dim);
                elseif sum_value < dim-1
                    % flipped around dim/2 face for rejection sampling
                    X = 1-DBMOPP.recalibrate(X,number_of_points,S,dim-sum_value,dim);
                else
                    % special case when sum_value >= dim-1, can just
                    % flip round the scaled unit simplex -- no rejection sampling
                    % needed :)
                    X = (X./repmat(S,1,dim))*(dim-sum_value);
                    X = 1-X;
                end
            end
            
        end
        
        function X = recalibrate(Z,number_of_points,S,sum_value,dim)
            X = (Z./repmat(S,1,dim))*sum_value;
            
            for i=1:number_of_points
                while (max(X(i,:)) > 1) % rejection sampling -- keep sampling until the constraints satisfied
                    Z(i,:) = exprnd(ones(1,dim));
                    S(i) = sum(Z(i,:));
                    X(i,:) = Z(i,:)/S(i) * sum_value;
                end
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
                    if (obj.inConvexHullOfAttractorRegion(x))
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
                    if (obj.inConvexHullOfAttractorRegion(x))
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
            if obj.getSoftConstraintViolation(x) > 0
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
                if obj.numberOfObjectives > 2
                    obj.attractorRegions{i}.convhull = convhull(locations(:,1),locations(:,2));
                end
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
                if (length(locations(:,1)) > 2) %only calculate the convex hull if there are more than two points
                    obj.attractorRegions{i}.convhull = convhull(locations(:,1),locations(:,2));
                end
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
            % set the neutral value to be the same in all neutral locations
            obj.neutralRegionObjectiveValues = ones(1,obj.numberOfObjectives)*obj.neutralRegionObjectiveValues;
            obj.placeVertexConstraintLocations();
            obj.placeCentreConstraintLocations();
            obj.placeMoatConstraintLocations();
            obj.assignDesignDimensionProjection();
            fprintf('set projection vectors\n');
            fprintf('set rescaling\n');
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
            % angles at offset1 and offset2 bracket that at the pivot
            
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
            inHardConstraintRegion = DBMOPP.inRegionExcludingBoundary(obj.hardConstraintCentres, obj.hardConstraintRadii, x);
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
                k = sum(d < obj.softConstraintRadii);
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
        function assignDesignDimensionProjection(obj)
            % if more than two design dimensions in problem, need to assign
            % the mapping down from this higher space to the 2D version
            % which will be subsequantly evaluated
            if obj.numberOfDesignVariables > 2
                if obj.variableSolutionDensity
                    difference = randi(obj.numberOfDesignVariables - 1);
                    mask = randperm(obj.numberOfDesignVariables);
                    mask = mask(1:difference);
                else
                    mask = randperm(obj.numberOfDesignVariables);
                    mask = mask(1:ceil(obj.numberOfDesignVariables/2));  % select random half of dimension indices
                end
                obj.pi1 = false(1,obj.numberOfDesignVariables);
                obj.pi1(1,mask) = true;
                obj.pi1Magnitude = sum(obj.pi1); % effectively norm(obj.pi1)^2; -- but quicker, and returns integer needed
                obj.pi2 = true(1,obj.numberOfDesignVariables);
                obj.pi2(1,mask) = false;
                obj.pi2Magnitude = sum(obj.pi2);
            end
        end
        %--
        function x = get2DVersion(obj,x)
            if length(x) > 2
                x = [dot(x,obj.pi1)/obj.pi1Magnitude, dot(x,obj.pi2)/obj.pi2Magnitude];
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
            I = find(d <= obj.centreRadii+eps);
            % check if inside or on any of the circles defining arractor regions
            index = -1;
            if isempty(I)==false
                % check if in Global Pareto region rather than e.g. local set
                if (I(1) > obj.numberOfLocalParetoSets) && (I(1) <= obj.numberOfLocalParetoSets + obj.numberOfGlobalParetoSets)
                    if (obj.constraintType == 2) || (obj.constraintType == 6)
                        % special case where Pareto set lies on the
                        % coundary of the circle due to the constraint
                        % covering the circle
                        if (abs(d(I(1))-obj.centreRadii(I(1))) < 1e4*eps(min(abs(d(I(1))),abs(obj.centreRadii(I(1)))))) % deal with rounding error
                            inHull = true;
                        end
                    elseif inpolygon(x(1),x(2),obj.attractorRegions{I(1)}.locations(obj.attractorRegions{I(1)}.convhull,1), obj.attractorRegions{I(1)}.locations(obj.attractorRegions{I(1)}.convhull,2))
                        % 'inside radius of global attractor region, need to check if in conv hull'
                        inHull = true;
                    end
                end
            end
            
            if obj.globalParetoSetType == 0 % identical performance
                inParetoRegion = inHull;
                inHull = false;
            elseif (obj.constraintType == 2) || (obj.constraintType == 6)
                inParetoRegion = inHull;
                inHull = false;
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

