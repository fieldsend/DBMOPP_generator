function [y] = distance_points_problem(x,distance_problem_parameters)
    
    if (length(x)>2)
       % projection being used, need to map back 
        x = project_nD_point_to_2D(x,distance_problem_parameters.projection_vectors(1,:),...
            distance_problem_parameters.projection_vectors(2,:));
    end
    
    y = zeros(1,length(distance_problem_parameters.num_objectives));
    for i=1:distance_problem_parameters.num_objectives
        d = minkowski_dist(x,distance_problem_parameters.distance_vectors(i).coordinates,distance_problem_parameters.minkowski_powers(i));
        d = modify_due_to_curvatue(d,x,distance_problem_parameters);
        
        % apply penalties
        for k=1:length(distance_problem_parameters.penalty_radii)
            pen_d = minkowski_dist(x,distance_problem_parameters.penalty_centre_list(k,:),2);
            if (pen_d < distance_problem_parameters.penalty_radii(k))
                d=d + distance_problem_parameters.penalty_radii(k);
            end
        end
        
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

function d = modify_due_to_curvatue(d,x,distance_problem_parameters)
    if (distance_problem_parameters.curvature_radius > 0)
        dc = minkowski_dist(x,distance_problem_parameters.centre_list,2);
        % get index of centre it is closed to, if any
        I = find(dc < distance_problem_parameters.radii*distance_problem_parameters.curvature_radius);
        if (length(I)>1)
            error('should not fall in more than one region');
        end
        if (length(I)==1)
        % increase value for any inside repulsion region
            d = d + min(distance_problem_parameters.radii)*distance_problem_parameters.curvature_radius;
        end
    end
end