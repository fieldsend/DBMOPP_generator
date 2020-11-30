function [X,Q] = gecco_2019_quiver(dpp,n)

% [X,Q] = gecco_2019_quiver(dpp,n)
%
% helper function for quiver plots
%
% dpp = problem instance
% n = grid resolution
%
% Jonathan Fieldsend, University of Exeter, 2018,2019
% See license information in package, available at 
% https://github.com/fieldsend/DBMOPP_generator


x = linspace(-1,1,n);
step = (x(2)-x(1))*0.8;
X = zeros(n,n);
k=1;
Q = zeros(4*n^2,2);
U = zeros(4*n^2,2);
counter = 1;
R = zeros(n,n,dpp.num_objectives);

% get qulities
for i=1:n
    for j=1:n
        t = distance_points_problem([x(i) x(j)],dpp);
        R(i,j,:) = t';
    end
end

for i=1:n
    for j=1:n
        number_domed = 0;
        if (i>1)
            if dominates(R(i-1,j,:),R(i,j,:))
               Q(counter,:) = [x(i) x(j)]';
               U(counter,:) = [-step 0]';
               counter = counter +1;
               number_domed = number_domed +1;
            end
            if (j>1)
                if dominates(R(i-1,j-1,:),R(i,j,:))
                    Q(counter,:) = [x(i) x(j)]';
                    U(counter,:) = [-step -step]';
                    counter = counter +1;
                    number_domed = number_domed +1;
                end
            end
            if (j<n)
                if dominates(R(i-1,j+1,:),R(i,j,:))
                    Q(counter,:) = [x(i) x(j)]';
                    U(counter,:) = [-step step]';
                    counter = counter +1;
                    number_domed = number_domed +1;
                end
            end
        end
        if (i<n)
            if dominates(R(i+1,j,:),R(i,j,:))
               Q(counter,:) = [x(i) x(j)]';
               U(counter,:) = [step 0]';
               counter = counter +1;
               number_domed = number_domed +1;
            end
            if (j>1)
                if dominates(R(i+1,j-1,:),R(i,j,:))
                    Q(counter,:) = [x(i) x(j)]';
                    U(counter,:) = [step -step]';
                    counter = counter +1;
                    number_domed = number_domed +1;
                end
            end
            if (j<n)
                if dominates(R(i+1,j+1,:),R(i,j,:))
                    Q(counter,:) = [x(i) x(j)]';
                    U(counter,:) = [step step]';
                    counter = counter +1;
                    number_domed = number_domed +1;
                end
            end
        end
        if (j>1)
            if dominates(R(i,j-1,:),R(i,j,:))
               Q(counter,:) = [x(i) x(j)]';
               U(counter,:) = [0 -step]';
               counter = counter +1;
               number_domed = number_domed +1;
            end
        end
        if (j<n)
            if dominates(R(i,j+1,:),R(i,j,:))
               Q(counter,:) = [x(i) x(j)]';
               U(counter,:) = [0 step]';
               counter = counter +1;
               number_domed = number_domed +1;
            end
        end
        if (number_domed > 0)
            % rescale
            U(counter-number_domed:counter-1,:)= U(counter-number_domed:counter-1,:)*number_domed;
            
        end
    end
    
    
end
Q(counter:end,:) = [];
    U(counter:end,:) = [];
    figure; 
    quiver(Q(:,1),Q(:,2),U(:,1),U(:,2));
    axis square;
    axis([-1 1 -1 1]);
end

function x = dominates(a,b)
    x = false;
    if (sum(a<b)>0) && (sum(a<=b)==length(a))
        x = true;
    end
end

        
