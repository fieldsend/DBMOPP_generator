function d = minkowski_dist(x,X,p)

% Function for arbitary Minkowski distances (NOT USED IN CURRENT GENERATOR)
%
% Jonathan Fieldsend, University of Exeter, 2018,2019
% See license information in package, available at 
% https://github.com/fieldsend/DBMOPP_generator


[n,m] = size(X);
[n1,m1] = size(x);
if (m~=m1)
   error('Arguments must have same column number'); 
end
if (n1~=1)
   error('First argument, x, must be a single row'); 
end

d = sum(abs(repmat(x,n,1)-X).^p,2).^(1/p);

end
