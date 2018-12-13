
function a = project_nD_point_to_2D(b,v1,v2)
a = [dot(b,v1)/(norm(v1)^2), dot(b,v2)/(norm(v2)^2)];

end