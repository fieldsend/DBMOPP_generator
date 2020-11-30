function a = project_nD_point_to_2D(b,v1,v2)

% helper function for converting between nD projection and 2D
%
% Jonathan Fieldsend, University of Exeter, 2018,2019
% See license information in package, available at 
% https://github.com/fieldsend/DBMOPP_generator


a = [dot(b,v1)/(norm(v1)^2), dot(b,v2)/(norm(v2)^2)];

end
