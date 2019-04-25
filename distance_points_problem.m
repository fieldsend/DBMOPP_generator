function [y] = distance_points_problem(x,distance_problem_parameters)
% function [y] = distance_points_problem(x,distance_problem_parameters)
%
% INPUTS
%
% x = design location, |x|>=2, must match design size expected in problem
%         instance
% distance_problem_parameters = problem instance (output from 
%         distance_point_generator function) [OPTIONAL, if not supplied
%         then a global variable named DISTANCE_PROBLEM_PARAMETERS must be 
%         set with this structure 
% OUTPUTS
%
% y = objective values assocaited with this x, given the problem instance
%
% Cost funtion takes in design vector and evaluates it under the problem 
% structure 
%
% Jonathan Fieldsend, University of Exeter, 2018,2019
% See license information in package, available at 
% https://github.com/fieldsend/DBMOPP_generator

global DISTANCE_PROBLEM_PARAMETERS

if (exist('distance_problem_parameters','var')==0)
    distance_problem_parameters = DISTANCE_PROBLEM_PARAMETERS;
end
if (length(x)>2)
    % projection being used, need to map back
    x = project_nD_point_to_2D(x,distance_problem_parameters.projection_vectors(1,:),...
        distance_problem_parameters.projection_vectors(2,:));
end

y = zeros(1,length(distance_problem_parameters.num_objectives));
for i=1:distance_problem_parameters.num_objectives
    % get distance to each objective attractor
    d = minkowski_dist(x,distance_problem_parameters.distance_vectors(i).coordinates,2);
    
    % apply any region penalties, which makes different pareto set regions non-identical -- potentially stacking
    for k=1:length(distance_problem_parameters.region_penalty_radii)
        pen_d = minkowski_dist(x,distance_problem_parameters.region_penalty_locations(k,:),2);
        if (pen_d < distance_problem_parameters.region_penalty_radii(k))
            d=d + distance_problem_parameters.region_penalty_radii(k);
        end
    end
    % apply any discontinuity and/or neutrality penalties
    
    % shift objective values
    d = d*distance_problem_parameters.objective_multiplier(i);
    d = d - distance_problem_parameters.objective_min(i);
    
    y(i) = min(d);
end

end



