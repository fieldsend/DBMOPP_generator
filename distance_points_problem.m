function [y] = distance_points_problem(x,distance_problem_parameters)
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
        d = minkowski_dist(x,distance_problem_parameters.distance_vectors(i).coordinates,2);
        
        % apply region penalties -- potentially stacking
        for k=1:length(distance_problem_parameters.region_penalty_radii)
            pen_d = minkowski_dist(x,distance_problem_parameters.region_penalty_locations(k,:),2);
            if (pen_d < distance_problem_parameters.region_penalty_radii(k))
                d=d + distance_problem_parameters.region_penalty_radii(k);
            end
        end
        % apply discontinuity penalties
        
        
        % apply neutrality penalties
        
        
        y(i) = min(d);
    end
    
end

function a = project_nD_point_to_2D(b,v1,v2)
a = [dot(b,v1)/(norm(v1)^2), dot(b,v2)/(norm(v2)^2)];

end

function d = minkowski_dist(x,X,p)

[n,~] = size(X);

d = sum(abs(repmat(x,n,1)-X).^p,2).^(1/p);

end
