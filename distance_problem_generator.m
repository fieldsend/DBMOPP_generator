function [distance_problem_parameters] = distance_problem_generator(num_objectives,num_dimensions,...
    curvature, number_of_disconnected_set_regions,...
    number_of_local_fronts,number_of_dominance_resistance_regions, ...
    number_of_discontinuous_regions,...
    varying_density,non_identical_pareto_sets,varying_objective_ranges, ...
    fill_space,plot_wanted,random_seed)

%[P] = distance_problem_generator(num_objectives,num_dimensions,curvature,...
%    number_of_disconnected_set_regions,number_of_local_fronts,...
%    number_of_dominance_reistance_regions,varying_objective_difficulty, ...
%    disconnect_pareto_front)
%
% num_objectives = number of objectives
% num_dimensions = number of dimensions [minium 2, OPTIONAL - default 2]
% curvature = flag to add concave region. [OPTIONAL, default false]
% number_of_disconnected_set_regions = number of disconnected set regions.
%                  [minium 1 (i.e. connected), OPTIONAL - default 1]
% number_of_local_fronts = number of additional local fronts. [minium 0 (i.e. single
%                  global), OPTIONAL - default 0]
% number_of_dominance_resistance_regions = number of dominance resistance
%                  regions [minium 0 (i.e. none), OPTIONAL - default 0]
% number_of_discontinuous_regions = number regions of discontinuarity (inc. 
%                  neutrality) [OPTIONAL - default 0]
% varying_density = flag to indicate if their can be a varying mapping
%                  density to x1 than x2 when num_dimensions>2 [OPTIONAL, 
%                  default false, must be false when num_dimensions=2] 
% non_identical_pareto_sets = flag to indicate if disconnected Pareto sets 
%                  are not identical (even when rotated back) [OPTIONAL,
%                  default false]
% varying_objective_ranges = flag to indicate objectives on different ranges
%                  [OPTIONAL - default false]
% fill_space = radii of regions set as large as possible [OPTIONAL -
%                  default false]
% plot_wanted = flag to indicate a plot of the problem is wanted [OPTIONAL -
%                  default false]
% random_seed = random seed to be used [OPTIONAL, current state of random
%                  generator used if not provided]
%
% Generator sets initial centre and angles to points via radius to place
% objective minimising coodinates. Uses an addition centre per extra region
% and rotation angle (3 parameters per addition region) on top of the 2+1+D
% for the Pareto set
%
% In the optional plot, red circles indicate circumference of zones for
% Pareto sets, black circles denote circumference of zones for local
% fronts, blue circles denote circumference of zones for dominance
% resistance regions. All points in green circles have penalties applied,
% causing discontinuities in objective landscaes, and cutting out regions
% of disconnected Pareto sets
%
% Jonathan Fieldsend, University of Exeter, 2018

% set up optional values
if (exist('num_dimensions','var')==0) || (num_dimensions < 2)
    num_dimensions = 2;
end

if (num_dimensions > 2)
    error('this implementation does not currently support design spaces > 2');
end

if (exist('curvature','var')==0)
    curvature = false;
end

if (curvature)
   error('Varying curvature not in this version \n'); 
end

if (exist('number_of_disconnected_set_regions','var')==0) || (number_of_disconnected_set_regions < 1)
    fprintf('Default used: single connected Pareto set\n');
    number_of_disconnected_set_regions = 1;
end

if (exist('number_of_local_fronts','var')==0) || (number_of_local_fronts < 0)
    fprintf('Default used: no local fronts\n');
    number_of_local_fronts = 0;
end

if (exist('number_of_dominance_resistance_regions','var')==0) || (number_of_dominance_resistance_regions < 0)
    fprintf('Default used: dominance resistance region number set at 0\n');
    number_of_dominance_resistance_regions = 0;
end

if (exist('number_of_discontinuous_regions','var')==0) || (number_of_discontinuous_regions < 0)
    fprintf('Default used: discontinuous region number set at 0\n');
    number_of_discontinuous_regions = 0;
end


if (number_of_discontinuous_regions>0)
   error('Arbitary discontinuos regions not in this version\n'); 
end

if exist('varying_density','var')==0
    fprintf('Default used: varying density set at false\n');
    varying_objective_difficulty = false;
end

if (varying_density) && (num_dimensions==2)
   error('Cannot vary density if design space is 2'); 
end

if exist('non_identical_pareto_sets','var')==0
    fprintf('Default used: Pareto sets fixed as identical\n');
    non_identical_pareto_sets = false;
end

if exist('varying_objective_ranges','var')==0
    fprintf('Default used: ranges of objective same and unchanges\n');
    varying_objective_ranges = false;
end

if exist('fill_space','var')==0
    fprintf('Default used: fill space set at false\n');
    fill_space=false; % apply random seed
end


if exist('plot_wanted','var')==0
    fprintf('Default used: no plot at end\n');
    plot_wanted=false; % apply random seed
end


if exist('random_seed','var')==1
    rng(abs(random_seed)); % apply random seed
end

% NOW ALLOCATE INSTANCE

number_of_penalty_regions = number_of_discontinuous_regions; %initial penalty locations
% number of different distinct regions that will need to be fitted
number_of_centres = number_of_disconnected_set_regions+number_of_local_fronts+...
    number_of_dominance_resistance_regions;

% set up coordinate vector holders
for i=1:num_objectives
    distance_problem_parameters.distance_vectors(i).coordinates=[];
end

% set up other holders
distance_problem_parameters.projection_vectors = []; % project from ND to 2D
distance_problem_parameters.num_objectives = num_objectives;
distance_problem_parameters.radii = zeros(number_of_centres,1); % radii of regions
distance_problem_parameters.curvature_radius = 0;
distance_problem_parameters.penalty_radii = zeros(number_of_discontinuous_regions,1);
distance_problem_parameters.objective_min = zeros(num_objectives,1);
distance_problem_parameters.objective_multiplier = ones(num_objectives,1);

if (varying_objective_ranges)
    % objective minimas between -100 and 100
    distance_problem_parameters.objective_min = (rand(num_objectives,1)-0.5)*200;
    % objective range multipliers between 1 and 1000
    distance_problem_parameters.objective_multiplier = rand(num_objectives,1)*1000;
end

% if want areas of concavity on Pareto front, use an internal circle with
% penalty. radius at random
if (curvature)
    distance_problem_parameters.curvature_radius = rand();
end


% arbitrary radius, scaled by the number of different distinct regions that
% will need to be fitted into the region
%
% ceil(1+2(N^0.5)) *2*r <= 2 [as between -1 and +1]

% set up parameters for Pareto regions, local fronts and dominance
% resistance regions
%radius  = 1.0/ceil(2*sqrt(number_of_centres + (number_of_discontinuous_regions/4.0))+1);
radius = 1.0/(ceil(sqrt(number_of_centres + number_of_discontinuous_regions/4))+2);
if (fill_space==false) % if not maximally filling space with regions
    radius = radius*rand(1,1); 
else
    fprintf('Radius set at max to fill space\n');
end
pareto_angles = get_random_angles(num_objectives); % arbitrary angles for Pareto set
distance_problem_parameters.rotations = rand(number_of_centres,1)*2*pi; %arbitary rotations for regions

%place centres of all regions (attractor and penalty)
[distance_problem_parameters,radius,radius_original] = place_regions(distance_problem_parameters,number_of_centres, number_of_penalty_regions, radius);
distance_problem_parameters.base_radius = radius;
fprintf('Centres assigned\n');
% assign local regions
centre_index=1;
[distance_problem_parameters, centre_index, radius] = ...
    assign_local_fronts(distance_problem_parameters,centre_index,number_of_local_fronts,pareto_angles,radius);
[distance_problem_parameters, centre_index] = ...
    assign_disconnected_set_regions(distance_problem_parameters,centre_index,number_of_disconnected_set_regions,pareto_angles,radius);
[distance_problem_parameters,region_types] = ...
    assign_penalty_regions_for_non_identical_sets(distance_problem_parameters,number_of_local_fronts,number_of_disconnected_set_regions,non_identical_pareto_sets);
[distance_problem_parameters, centre_index] = ...
    assign_dominance_resistance_regions(distance_problem_parameters,centre_index,number_of_dominance_resistance_regions,pareto_angles,radius,region_types);
[distance_problem_parameters] = ...
    assign_penalty_regions(distance_problem_parameters,centre_index,number_of_discontinuous_regions,radius_original);

if (plot_wanted)
    % now plot problem instance
    plot_2D_regions(distance_problem_parameters,number_of_centres,number_of_discontinuous_regions,number_of_local_fronts,number_of_disconnected_set_regions,radius,num_objectives);
end
end


% ---- HELPER FUNCTIONS ---

function plot_2D_regions(distance_problem_parameters,...
    number_of_centres,number_of_discontinuous_regions,number_of_local_fronts,number_of_disconnected_set_regions,radius,num_objectives)

close all; figure; hold on;
plot(distance_problem_parameters.centre_list(:,1), distance_problem_parameters.centre_list(:,2),'kx');
%plot(centre_list(number_of_local_fronts:number_of_local_fronts+number_of_disconnected_set_regions-1,1), centre_list(number_of_local_fronts:number_of_local_fronts+number_of_disconnected_set_regions-1,2),'r*');
axis([-1 1 -1 1])
for i=1:number_of_centres-number_of_discontinuous_regions
    if i<=number_of_local_fronts
        c='k-';
        add_circle_to_plot(distance_problem_parameters.centre_list(i,:),radius,c);
    elseif i <= number_of_local_fronts+number_of_disconnected_set_regions
        c='r-';
        add_circle_to_plot(distance_problem_parameters.centre_list(i,:),radius,c);
    else % dominance resistance regions
        c='b-';
        add_circle_to_plot(distance_problem_parameters.centre_list(i,:),radius,c);
        %add_circle_to_plot(centre_list(i,:),4*radius,c);
    end
    
end
for i=1:number_of_discontinuous_regions
    add_circle_to_plot(distance_problem_parameters.penalty_centre_list(i,:),distance_problem_parameters.penalty_radii(i),'g-');
end


for i=1:length(distance_problem_parameters.region_penalty_radii)
    add_circle_to_plot(distance_problem_parameters.region_penalty_locations(i,:),distance_problem_parameters.region_penalty_radii(i),'g-');
end

axis square

for i=1:num_objectives
    plot(distance_problem_parameters.distance_vectors(i).coordinates(:,1),distance_problem_parameters.distance_vectors(i).coordinates(:,2),'bo');
end

xlabel('x_1')
ylabel('x_2')

end


function angles = get_random_angles(n)
angles = rand(n,1)*2*pi; % arbitrary angles for points on circulference of region
end

function [distance_problem_parameters, centre_index, radius] = assign_local_fronts(distance_problem_parameters,centre_index,number_of_local_fronts,angles,radius)
if (number_of_local_fronts > 0)
    % need to rescale radius as local fronts need to have wider sets
    radius = radius/2;
    scaling = linspace(1,2,number_of_local_fronts+1);
    for j=2:number_of_local_fronts+1
        V = repmat(distance_problem_parameters.centre_list(centre_index,:),length(angles),1)...
            +repmat(radius*scaling(j),1,2).*[cos(angles+distance_problem_parameters.rotations(centre_index)), sin(angles+distance_problem_parameters.rotations(centre_index))];
        for k=1:length(angles)
            distance_problem_parameters.distance_vectors(k).coordinates = ...
                [distance_problem_parameters.distance_vectors(k).coordinates; V(k,:)];
        end
        distance_problem_parameters.radii(centre_index) = radius*scaling(j);
        centre_index = centre_index+1;
    end
end
end

function [distance_problem_parameters, centre_index] = assign_disconnected_set_regions(distance_problem_parameters,centre_index,number_of_disconnected_set_regions,angles,radius)
% assign disconnected sets
for i=1:number_of_disconnected_set_regions
    V = repmat(distance_problem_parameters.centre_list(centre_index,:),length(angles),1)+...
        repmat(radius,1,2).*[cos(angles+distance_problem_parameters.rotations(centre_index)), sin(angles+distance_problem_parameters.rotations(centre_index))];
    for k=1:length(angles)
        distance_problem_parameters.distance_vectors(k).coordinates = [distance_problem_parameters.distance_vectors(k).coordinates; V(k,:)];
    end
    distance_problem_parameters.radii(centre_index) = radius;
    centre_index = centre_index+1;
end
end

function [distance_problem_parameters, centre_index] = assign_dominance_resistance_regions(distance_problem_parameters,centre_index,number_of_dominance_resistance_regions,angles,radius,region_types)
n = length(angles);
% assign dominance resistance regions
for i=1:number_of_dominance_resistance_regions
    V = repmat(distance_problem_parameters.centre_list(centre_index,:),length(angles),1)...
        +repmat(radius,1,2).*[cos(angles+distance_problem_parameters.rotations(centre_index)), sin(angles+distance_problem_parameters.rotations(centre_index))];
    if isempty(region_types)==false
        % special case where non-identical pareto sets used, so penalties
        % need to also be applied to dominance resistance locations
        [region_penalty_locations,region_penalty_radii] = apply_local_front_penalties_dr(V,region_types);
        distance_problem_parameters.region_penalty_locations = [distance_problem_parameters.region_penalty_locations; region_penalty_locations];
        distance_problem_parameters.region_penalty_radii = [distance_problem_parameters.region_penalty_radii; region_penalty_radii];
    end
    % now flag which objectives to include
    num_to_include = randperm(n-1);
    num_to_include = num_to_include(1); % number to include
    [~,I] = sort(rand(n,1)); %random objective indices
    for k=1:num_to_include
        distance_problem_parameters.distance_vectors(I(k)).coordinates = ...
            [distance_problem_parameters.distance_vectors(I(k)).coordinates; V(I(k),:)];
    end
    distance_problem_parameters.radii(centre_index) = radius;
    centre_index = centre_index+1;
end

end

function [distance_problem_parameters,region_types] = assign_penalty_regions_for_non_identical_sets(distance_problem_parameters,number_of_local_fronts,number_of_disconnected_set_regions,non_identical_pareto_sets)
distance_problem_parameters.set_penalty_types = 0;
distance_problem_parameters.region_penalty_locations = [];
distance_problem_parameters.region_penalty_radii = [];
region_types=[];
if (non_identical_pareto_sets)
    if (number_of_disconnected_set_regions>1)
        % if disconnected Pareto sets, and what to be non-identical
        %local fronts allocated first, then disconnected regions
        
        % as many distinct penalty types as there are disconnected regions
        distance_problem_parameters.set_penalty_types = number_of_disconnected_set_regions;
        % different subset of objectives of front penalised with each 
        
        region_index=number_of_local_fronts+1;
        region_types=[];
        for i=1:number_of_disconnected_set_regions
            [region_types{i}.region_penalty_locations,region_types{i}.region_penalty_radii, region_types{i}.affected_indices] =...
                randomly_apply_region_penalties(distance_problem_parameters,region_index);
            region_index = region_index+1;
            distance_problem_parameters.region_penalty_locations = [distance_problem_parameters.region_penalty_locations; region_types{i}.region_penalty_locations];
            distance_problem_parameters.region_penalty_radii = [distance_problem_parameters.region_penalty_radii; region_types{i}.region_penalty_radii];
        end
        % region types holds the distinct penalty forms used, which must be
        % selected from to apply to each local front to ensure none contact
        % pareto optimal solutions
        region_index=1;
        for i=1:number_of_local_fronts
            [region_penalty_locations,region_penalty_radii] = apply_local_front_penalties(distance_problem_parameters,region_index,region_types);
            region_index = region_index+1;
            distance_problem_parameters.region_penalty_locations = [distance_problem_parameters.region_penalty_locations; region_penalty_locations];
            distance_problem_parameters.region_penalty_radii = [distance_problem_parameters.region_penalty_radii; region_penalty_radii];
        end
    end
end

end

function [region_penalty_locations,region_penalty_radii,affected_indices] = randomly_apply_region_penalties(distance_problem_parameters,region_index)
    affected_indices = randperm(length(distance_problem_parameters.distance_vectors));
    max_I = randperm(length(distance_problem_parameters.distance_vectors)-1);
    max_I = max_I(1);
    affected_indices = affected_indices(1:max_I);
    % now selected 1 to number_objectives-1 locations to add penalties,
    % width of penaltie radius needs to be proportional to radius of base
    % region
    region_penalty_radii = distance_problem_parameters.base_radius*rand(length(affected_indices),1)/2;
    region_penalty_locations = zeros(length(affected_indices),2);
    for i=1:length(affected_indices)
        % append minima location of region
        region_penalty_locations(i,:) = distance_problem_parameters.distance_vectors(affected_indices(i)).coordinates(region_index,:);
    end
end        

function [region_penalty_locations,region_penalty_radii] = apply_local_front_penalties(distance_problem_parameters,region_index,region_types)
    I = randperm(length(region_types));
    I = I(1);
    region_penalty_locations = zeros(length(region_types{I}.affected_indices),2);
    region_penalty_radii = region_types{I}.region_penalty_radii*(distance_problem_parameters.radii(region_index)/distance_problem_parameters.base_radius); %rescale by local front multiplier
    
    for i=1:length(region_types{I}.affected_indices)
        % append minima location of region which mimic one of the sets of
        % penalties already applied to a Pareto set region
        region_penalty_locations(i,:) = distance_problem_parameters.distance_vectors(region_types{I}.affected_indices(i)).coordinates(region_index,:);
    end    

end

function [region_penalty_locations,region_penalty_radii] = apply_local_front_penalties_dr(V,region_types)
    I = randperm(length(region_types));
    I = I(1);
    region_penalty_locations = zeros(length(region_types{I}.affected_indices),2);
    region_penalty_radii = region_types{I}.region_penalty_radii;
    for i=1:length(region_types{I}.affected_indices)
        % append minima location of region which mimic one of the sets of
        % penalties already applied to a Pareto set region
        region_penalty_locations(i,:) = V(region_types{I}.affected_indices(i),:);
    end
end
                


function [distance_problem_parameters, centre_index] = assign_penalty_regions(distance_problem_parameters,centre_index,number_of_discontinuous_regions,radius_original)

% assign penalty regions
for i=1:number_of_discontinuous_regions
    distance_problem_parameters.penalty_radii(i) = rand(1,1)*radius_original;
end
end


function [distance_problem_parameters,radius,radius_original] = place_regions(distance_problem_parameters,number_of_centres, number_of_penalty_regions, radius)
% number of centres required for each region
distance_problem_parameters.centre_list = zeros(number_of_centres,2);
radius_original=radius; % keep track of max radius if local fronts used
% rejection sampling to place regions

% place pareto, local and dominance resiatant regions first, need to be
% 4r apart from each other centre
time_start = tic; 
time_elapsed = 0;
MAX_ELAPSED = 5; % max seconds before reattempting
invalid = true;
while (invalid)
    distance_problem_parameters.centre_list(1,:) = (rand(1,2)*2)-1; % arbitrary centre of Pareto set
    if ((sum((distance_problem_parameters.centre_list(1,:)+radius)>1)==0) && (sum((distance_problem_parameters.centre_list(1,:)-radius)<-1)==0))
        invalid = false; % ensure set doesn't cross boundary of feasible space
    end
end
fprintf('Radius: %f\n',radius);
% now assign centre locations
for i=2:number_of_centres
    invalid = true;
    while(invalid)
        r = (rand(1,2)*2)-1; %random cordinate pair between -1 and 1
        if ((sum((r+radius)>1)==0) && (sum((r-radius)<-1)==0)) %ensure not too close to boundary
            t = min(minkowski_dist(r, distance_problem_parameters.centre_list(1:i-1,:),2));
            if ( t > 4*radius) % centres not equal or closer than 4*radius
                fprintf('Assigned centre %d\n', i);
                invalid=false;
            end
        end
        time_elapsed = toc(time_start);
        if (time_elapsed>MAX_ELAPSED)
            break;
        end
    end
    distance_problem_parameters.centre_list(i,:) = r;
end
% now place penalty regions for discontinuous and neutral objective
% regions, need to be at least r from each other centre
distance_problem_parameters.penalty_centre_list = zeros(number_of_penalty_regions,2);
for i=2:number_of_penalty_regions
    invalid = true;
    while(invalid)
        r = (rand(1,2)*2)-1;
        if ((sum((r+radius)>1)==0) && (sum((r-radius)<-1)==0))
            t = min(minkowski_dist(r, [distance_problem_parameters.centre_list; distance_problem_parameters.penalty_centre_list(1:i-1,:)],2));
            if ( t > 2*radius) % centres not equal or closer than 2*radius
                fprintf('Assigned penalty centre %d\n', i);
                invalid=false;
            end
        end
        time_elapsed = toc(time_start);
        if (time_elapsed>MAX_ELAPSED)
            break;
        end
    end
    distance_problem_parameters.penalty_centre_list(i,:) = r;
end

if (time_elapsed>MAX_ELAPSED)
    fprintf('restarting...\n');
    [distance_problem_parameters,radius,radius_original] = place_regions(distance_problem_parameters,number_of_centres, number_of_penalty_regions, radius*0.95);
end


end

function add_circle_to_plot(centre,radius,c)
step = 0:pi/100:2*pi;
x = radius * cos(step) + centre(1);
y = radius * sin(step) + centre(2);
plot(x, y,c);
end