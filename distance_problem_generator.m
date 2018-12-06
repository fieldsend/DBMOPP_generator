function [distance_problem_parameters] = distance_problem_generator(num_objectives,num_dimensions,...
    curvature, number_of_disconnected_set_regions,...
    number_of_local_fronts,number_of_dominance_resistance_regions, ...
    number_of_discontinuous_regions,...
    varying_objective_difficulty,disconnected_pareto_front,random_seed)

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
% number_of_local_fronts = number of local fronts. [minium 1 (i.e. single 
%                  global), OPTIONAL - default 1]
% number_of_dominance_resistance_regions = number of dominance resistance 
%                  regions [minium 0 (i.e. none), OPTIONAL - default 0]
% varying_objective_difficulty = effective only for num_dimensions > 2 
%                   [OPTIONAL - default false]
% disconnected_pareto_front = use reuplsion points to create disconnected 
%                   Pareto fronts [OPTIONAL - default false]
%
% Generator sets initial centre and angles to points via radius to place 
% objective minimising coodinates. Uses an addition centre per extra region
% and rotation angle (3 parameters per addition region) on top of the 2+1+D
% for the Pareto set
%
% Jonathan Fieldsend, University of Exeter, 2018

if (exist('num_dimensions','var')==0) || (num_dimensions < 2)
    num_dimensions = 2;
end

if (exist('curvature','var')==0)
    curvature = false;
end


if (exist('number_of_disconnected_set_regions','var')==0) || (number_of_disconnected_set_regions < 1)
    number_of_disconnected_set_regions = 1;
end

if (exist('number_of_local_fronts','var')==0) || (number_of_local_fronts < 1)
    number_of_local_fronts = 1;
end

if (exist('number_of_dominance_resistance_regions','var')==0) || (number_of_dominance_resistance_regions < 0)
    number_of_dominance_resistance_regions = 0;
end

if (exist('number_of_discontinuous_regions','var')==0) || (number_of_discontinuous_regions < 0)
    number_of_discontinuous_regions = 0;
end

if exist('varying_objective_difficulty','var')==0
    varying_objective_difficulty = false;
end

if exist('disconnected_pareto_front','var')==0
    disconnected_pareto_front = false;
end

% number of different distinct regions that will need to be fitted
number_of_centres = number_of_disconnected_set_regions+number_of_local_fronts-1+...
    number_of_dominance_resistance_regions+number_of_discontinuous_regions;
for i=1:num_objectives
    distance_problem_parameters.distance_vectors(i).coordinates=[];
end
distance_problem_parameters.minkowski_powers = ones(num_objectives,1)*2;
distance_problem_parameters.repulsion_vectors = [];
distance_problem_parameters.projection_vectors = [];
distance_problem_parameters.num_objectives = num_objectives;
distance_problem_parameters.radii = zeros(number_of_centres-number_of_discontinuous_regions,1);
distance_problem_parameters.curvature_radius = 0;
distance_problem_parameters.penalty_radii = zeros(number_of_discontinuous_regions,1);


% if want areas of concavity on Pareto front, use an internal cicrle with
% penalty. radius at random
if(curvature)
    distance_problem_parameters.curvature_radius = rand();
end

% arbitrary radius, scaled by the number of different distinct regions that
% will need to be fitted into the region
%
% ceil(1+2(N^0.5)) *2*r <= 2

radius = (rand(1,1))/ceil(2*sqrt(number_of_centres)+1);
radius_m = radius;
%radius = 1/ceil(2*sqrt(number_of_centres)+1);
angles = rand(num_objectives,1)*2*pi; % arbitrary angles for Pareto set
rotations = rand(number_of_centres,1)*2*pi; %arbitary rotations for regions
min_seperation = radius*4; %threshold for min distance between region centres

% number of centres required for each region
centre_list = zeros(number_of_centres,2);

invalid = true;
while (invalid)
    centre_list(1,:) = (rand(1,2)*2)-1; % arbitrary centre of Pareto set
    if ((sum((centre_list(1,:)+radius)>1)==0) && (sum((centre_list(1,:)-radius)<-1)==0))
       invalid = false; % ensure set doesn't cross boundary of feasible space 
    end
end
fprintf('Radius: %f\n',radius);
% now assign centre locations
for i=2:number_of_centres
    invalid = true;
    while(invalid)
        r = (rand(1,2)*2)-1;
        if ((sum((r+radius)>1)==0) && (sum((r-radius)<-1)==0))
            t = min(sqrt(dist2(r, centre_list(1:i-1,:))));
            if ( t > min_seperation)
                fprintf('Assigned centre %d\n', i);
                invalid=false; 
            end
        end
    end
    centre_list(i,:) = r;
end
penalty_centre_list = centre_list(number_of_centres-number_of_discontinuous_regions+1:end,:);
centre_list = centre_list(1:number_of_centres-number_of_discontinuous_regions,:);

fprintf('Centres assigned\n');
% set minkowski power
% if vary_curvature
%    for i = 1:num_objectives
%        r = rand();
%        if r<1/3
%            distance_problem_parameters.minkowski_powers(i) = rand();
%        elseif r<2/3
%            distance_problem_parameters.minkowski_powers(i) = 1;
%        else
%            distance_problem_parameters.minkowski_powers(i) = 1+rand()*10;
%        end
%    end
% else
%     if curvature == -1
%         distance_problem_parameters.minkowski_powers = ones(num_objectives,1)*rand(); 
%     elseif curvature == 1
%         distance_problem_parameters.minkowski_powers = distance_problem_parameters.minkowski_powers + rand()*10;
%     end
% end

% assign local front regions
centre_index=1;
if (number_of_local_fronts > 1)
    % need to rescale radius as local fronts need to have wider sets
    radius = radius/2;
    scaling = linspace(1,2,number_of_local_fronts);
    for j=2:number_of_local_fronts
        V = repmat(centre_list(centre_index,:),length(angles),1)+repmat(radius*scaling(j),1,2).*[cos(angles+rotations(centre_index)), sin(angles+rotations(centre_index))];
        for k=1:num_objectives
            distance_problem_parameters.distance_vectors(k).coordinates = [distance_problem_parameters.distance_vectors(k).coordinates; V(k,:)];
        end
        distance_problem_parameters.radii(centre_index) = radius*scaling(j);
        centre_index = centre_index+1; 
    end
end
% assign disconnected sets
for i=1:number_of_disconnected_set_regions
    V = repmat(centre_list(centre_index,:),length(angles),1)+repmat(radius,1,2).*[cos(angles+rotations(centre_index)), sin(angles+rotations(centre_index))];
    for k=1:num_objectives
        distance_problem_parameters.distance_vectors(k).coordinates = [distance_problem_parameters.distance_vectors(k).coordinates; V(k,:)];
    end
    distance_problem_parameters.radii(centre_index) = radius;
    centre_index = centre_index+1; 
end
% assign dominance resistance regions
for i=1:number_of_dominance_resistance_regions
    V = repmat(centre_list(centre_index,:),length(angles),1)+repmat(radius,1,2).*[cos(angles+rotations(centre_index)), sin(angles+rotations(centre_index))];
    % now flag which objectives to remove
    num_to_exclude = randperm(num_objectives-1);
    num_to_exclude = num_to_exclude(1); % number to exclude
    [~,I] = sort(rand(num_objectives,1));
    for k=1:num_to_exclude
        distance_problem_parameters.distance_vectors(I(k)).coordinates = [distance_problem_parameters.distance_vectors(I(k)).coordinates; V(I(k),:)];
    end
    distance_problem_parameters.radii(centre_index) = radius;
    centre_index = centre_index+1; 
end
% assign penalty regions
for i=1:number_of_discontinuous_regions
   distance_problem_parameters.penalty_radii(i) = rand(1,1)*radius_m;
end



close all; figure; hold on;
plot(centre_list(:,1), centre_list(:,2),'kx');
%plot(centre_list(number_of_local_fronts:number_of_local_fronts+number_of_disconnected_set_regions-1,1), centre_list(number_of_local_fronts:number_of_local_fronts+number_of_disconnected_set_regions-1,2),'r*');
axis([-1 1 -1 1])
for i=1:number_of_centres-number_of_discontinuous_regions
    c='r-';
    if i<number_of_local_fronts
        c='k-';
        add_circle_to_plot(centre_list(i,:),radius,c);
    elseif i >= number_of_local_fronts+number_of_disconnected_set_regions
        c='k-';
        add_circle_to_plot(centre_list(i,:),radius,c);
    else
        add_circle_to_plot(centre_list(i,:),radius,c);
        %add_circle_to_plot(centre_list(i,:),4*radius,c);
    end
    
end
for i=1:number_of_discontinuous_regions
    add_circle_to_plot(penalty_centre_list(i,:),distance_problem_parameters.penalty_radii(i),'g-');
end

axis square

for i=1:num_objectives
    plot(distance_problem_parameters.distance_vectors(i).coordinates(:,1),distance_problem_parameters.distance_vectors(i).coordinates(:,2),'bo');
end

distance_problem_parameters.centre_list = centre_list; %return centre list
distance_problem_parameters.penalty_centre_list=penalty_centre_list;
xlabel('x_1')
ylabel('x_2')

% front curvature
% can vary minkowski power to alter from concave to convex between
% objectives


% disconnectness of Pareto set
% can place identical layouts with different centres

% disconnectness of Pareto front
% can place inversion points (objective *increases* as you get close)



% local fronts
% can place other centres with varying radii of Pareto set


% dominance resistance
% can place identical layouts with different centres and objective subsets

% correlation between objectives
% need to perform search on this....

% difficult/easy criteria
% vary projection vector length
end

function add_circle_to_plot(centre,radius,c)
    step = 0:pi/100:2*pi;
    x = radius * cos(step) + centre(1);
    y = radius * sin(step) + centre(2);
    plot(x, y,c);
end