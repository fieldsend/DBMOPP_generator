function plot_dbmopp_2D_regions(distance_problem_parameters,...
    number_of_discontinuous_regions,num_objectives,number_of_local_fronts,num_dom_res_regions)

% function plot_dbmopp_2D_regions(distance_problem_parameters,...
%    number_of_discontinuous_regions,num_objectives,number_of_local_fronts,num_dom_res_regions)
%
% Used by the generator to plot some of the problem attractors/penalty
% locations
%
% Jonathan Fieldsend, University of Exeter, 2018,2019
% See license information in package, available at 
% https://github.com/fieldsend/DBMOPP_generator
    

number_of_centres = size(distance_problem_parameters.centre_list,1);
radius = min(distance_problem_parameters.radii);



close all; figure; hold on;
%plot(distance_problem_parameters.centre_list(:,1), distance_problem_parameters.centre_list(:,2),'kx');
%plot(centre_list(number_of_local_fronts:number_of_local_fronts+number_of_disconnected_set_regions-1,1), centre_list(number_of_local_fronts:number_of_local_fronts+number_of_disconnected_set_regions-1,2),'r*');
axis([-1 1 -1 1])
for i=number_of_local_fronts+1:number_of_centres-num_dom_res_regions
    plot(distance_problem_parameters.centre_list(i,1), distance_problem_parameters.centre_list(i,2),'k+');
    add_circle_to_plot(distance_problem_parameters.centre_list(i,:),radius,'k-');
end
% for i=1:number_of_centres-number_of_discontinuous_regions
%     if i<=number_of_local_fronts
%         c='k-';
%         add_circle_to_plot(distance_problem_parameters.centre_list(i,:),radius,c);
%     elseif i <= number_of_local_fronts+number_of_disconnected_set_regions
%         c='r-';
%         add_circle_to_plot(distance_problem_parameters.centre_list(i,:),radius,c);
%     else % dominance resistance regions
%         c='b-';
%         add_circle_to_plot(distance_problem_parameters.centre_list(i,:),radius,c);
%         %add_circle_to_plot(centre_list(i,:),4*radius,c);
%     end
%     
% end
for i=1:number_of_discontinuous_regions
    add_circle_to_plot(distance_problem_parameters.penalty_centre_list(i,:),distance_problem_parameters.penalty_radii(i),'r-');
end


for i=1:length(distance_problem_parameters.region_penalty_radii)
    add_circle_to_plot(distance_problem_parameters.region_penalty_locations(i,:),distance_problem_parameters.region_penalty_radii(i),'r-');
end

axis square
my_col = jet(num_objectives);
for i=1:num_objectives
    plot(distance_problem_parameters.distance_vectors(i).coordinates(:,1),distance_problem_parameters.distance_vectors(i).coordinates(:,2),'.','color',my_col(i,:),'MarkerSize',20);
end

xlabel('x_1')
ylabel('x_2')

end



function add_circle_to_plot(centre,radius,c)
step = 0:pi/100:2*pi;
x = radius * cos(step) + centre(1);
y = radius * sin(step) + centre(2);
plot(x, y,c);
end
