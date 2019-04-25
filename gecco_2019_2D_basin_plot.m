function [perimeter_list, optima_list, region_list,mode_matrix,basin_matrix, basin_boundary_matrix] ...
    = gecco_2019_2D_basin_plot(dpp,n)
 
% function [perimeter_list, optima_list, region_list,mode_matrix,basin_matrix, basin_boundary_matrix] ...
%    = gecco_2019_2D_basin_plot(dpp,n)   
%
%dpp = distance-based point structure from generator
% n = number of samples per axis (grid resolution)
%
% Generates various 2D basin plots for the given problem structue at the 
% grid resolution provided
%
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

% get qualities for each
for i=1:n
    for j=1:n
        t = distance_points_problem([x(i), x(j)],dpp);
        R(i,j,:) = t';
    end
end

% now identify local optima regions
[perimeter_list, optima_list,region_list,mode_matrix,basin_matrix, basin_boundary_matrix] = get_perimeter_and_optima(R,n);



figure; plot((perimeter_list(:,1)/(n/2))-1,(perimeter_list(:,2)/(n/2))-1,'k.');
axis([-1 1 -1 1;]);
axis square;
figure; plot((optima_list(:,1)/(n/2))-1,(optima_list(:,2)/(n/2))-1,'r.');
axis([-1 1 -1 1;]);
axis square;

figure; imagesc(mode_matrix')
set(gca,'YDir','normal')

figure; imagesc(basin_matrix')
set(gca,'YDir','normal')


figure; imagesc(basin_boundary_matrix')
set(gca,'YDir','normal')

B = basin_boundary_matrix;
for i=1:n
    for j=1:n
        if (mode_matrix(i,j)~=0)
            B(i,j)=-1;
        end
    end
end
figure; imagesc(B')
set(gca,'YDir','normal')

% figure;
% hold on;
% for i=1:length(perimeter_list)
%     text(perimeter_list(i,1),perimeter_list(i,2), int2str(region_list(i)));
% end


end

function x = dominates(a,b)
x = false;
if (sum(a<b)>0) && (sum(a<=b)==length(a))
    x = true;
end
end


function [perimeter_list, optima_list, region_list,mode_matrix,basin_matrix,basin_boundary_matrix] = get_perimeter_and_optima(R,n)

perimeter_list = zeros(0,2);
optima_list = zeros(0,2);
mode_matrix = zeros(n,n);
basin_matrix = zeros(n,n);
perimeter_matrix = zeros(n,n);
basin_boundary_matrix = zeros(n,n);
% find modal members and fill optima_list
% using moore neighbourhood
for i=1:n
    for j=1:n
        number_dominated_by = 0;
        if (i>1)
            if dominates(R(i-1,j,:),R(i,j,:))
                number_dominated_by = number_dominated_by +1;
            end
            if (j>1)
                if dominates(R(i-1,j-1,:),R(i,j,:))
                    number_dominated_by = number_dominated_by +1;
                end
            end
            if (j<n)
                if dominates(R(i-1,j+1,:),R(i,j,:))
                    number_dominated_by = number_dominated_by +1;
                end
            end
        end
        if (i<n)
            if dominates(R(i+1,j,:),R(i,j,:))
                number_dominated_by = number_dominated_by +1;
            end
            if (j>1)
                if dominates(R(i+1,j-1,:),R(i,j,:))
                    number_dominated_by = number_dominated_by +1;
                end
            end
            if (j<n)
                if dominates(R(i+1,j+1,:),R(i,j,:))
                    number_dominated_by = number_dominated_by +1;
                end
            end
        end
        if (j>1)
            if dominates(R(i,j-1,:),R(i,j,:))
                number_dominated_by = number_dominated_by +1;
            end
        end
        if (j<n)
            if dominates(R(i,j+1,:),R(i,j,:))
                number_dominated_by = number_dominated_by +1;
            end
        end
        if (number_dominated_by == 0)
            optima_list(end+1,:) = [i,j]; % no neighbours dominate location i,j
            mode_matrix(i,j)=1;
        end
    end
end
% determine neighbours of each cell
Neighbours=cell(n,n);
dominating_neighbours=cell(n,n);
for i=1:n
    for j=1:n
        indices =[];
        if (i>1)
            indices(end+1,:) = [i-1, j];
            if (j>1)
                indices(end+1,:) = [i-1, j-1];
            end
            if (j<n)
                indices(end+1,:) = [i-1, j+1];
            end
        end
        if (i<n)
            indices(end+1,:) = [i+1, j];
            if (j>1)
                indices(end+1,:) = [i+1, j-1];
            end
            if (j<n)
                indices(end+1,:) = [i+1, j+1];
            end
        end
        if (j>1)
            indices(end+1,:) = [i, j-1];
        end
        if (j<n)
            indices(end+1,:) = [i, j+1];
        end
        Neighbours{i,j}=indices;
        dominating_neighbours{i,j}=[];
        for k=1:length(Neighbours{i,j})
            if (dominates(R(indices(k,1),indices(k,2),:),R(i,j,:)))
                dominating_neighbours{i,j}=[dominating_neighbours{i,j}; [Neighbours{i,j}(k,1) Neighbours{i,j}(k,2)]];
            end
        end
    end
end

%now determine perimeter list
for k=1:length(optima_list)
    i = optima_list(k,1);
    j = optima_list(k,2);
    perimeter_flag = true;
    for m=1:length(Neighbours{i,j})
        if (mode_matrix(Neighbours{i,j}(m,1),Neighbours{i,j}(m,2))==0)
            perimeter_list(end+1,:) = [i j];
            mode_matrix(i,j) = 2;
            break;
        end
    end
end

fprintf('iterate\n');
% iterative determine region -- worth checking if quicker way!

region_list = 1:length(perimeter_list)';
changed = true;
D = dist2(perimeter_list,perimeter_list);
D1 = dist2(perimeter_list(:,1),perimeter_list(:,1));
D2 = dist2(perimeter_list(:,2),perimeter_list(:,2));
D(D==0) = inf; % set diagonal to infinity
step = 0;
while(changed==true)
    [region_list,changed] = merge_regions(perimeter_list,region_list,D,D1,D2);
    step=step+1;
    fprintf('iteration %d, number of regions %d\n',step, length(unique(region_list)));
end

% label perimeter
for k=1:length(perimeter_list)
    perimeter_matrix(perimeter_list(k,1),perimeter_list(k,2)) = region_list(k);
end
basin_matrix_list = cell(n,n);    
for i=1:n
    fprintf('matrix row %d processing... \n',i);
    for j=1:n
        if (mode_matrix(i,j)==0) % if not a mode
            [attractor_list, basin_matrix_list] = recursive_attractor([i j], dominating_neighbours, mode_matrix,perimeter_matrix,basin_matrix_list);
            a = unique(attractor_list);
            basin_matrix(i,j) = length(a);
        end
    end
end
for i=1:n
    for j=1:n
        a = basin_matrix_list{i,j};
        if ((mode_matrix(i,j)==0) || (length(a)>1)) % if not a mode, or if attracted to multiple optima
            same = true;
            for k=1:length(Neighbours{i,j})
                a1 = basin_matrix_list{Neighbours{i,j}(k,1),Neighbours{i,j}(k,2)};
                if (isempty(a1))==true % don't draw boundary around neutral region
                    if perimeter_matrix(Neighbours{i,j}(k,1),Neighbours{i,j}(k,2))~=a % if region on perimeter isn't the attractor
                        same = false;
                    end
                else
                    if(length(a1)~=1)
                        same = false;
                    elseif (a~=a1)
                        same = false;
                    end
                end
            end
            if (same==false) 
                basin_boundary_matrix(i,j)=1;
            end
        end
    end
end


end

function [attractor_list,basin_matrix_list]  = recursive_attractor(x, dominating_neighbours, mode_matrix,perimeter_matrix,basin_matrix_list)
    attractor_list = [];
    [n,~] = size(dominating_neighbours{x(1),x(2)});
    if (mode_matrix(x(1),x(2))~=0) % not located at mode
        attractor_list =  perimeter_matrix(x(1),x(2));
    elseif isempty(basin_matrix_list{x(1),x(2)})==false
        attractor_list = basin_matrix_list{x(1),x(2)}; % precomputed, so use!
    else
        for i=1:n
             [a,basin_matrix_list] = recursive_attractor(dominating_neighbours{x(1),x(2)}(i,:), dominating_neighbours, mode_matrix,perimeter_matrix,basin_matrix_list);
             attractor_list = [attractor_list a];
        end
        attractor_list = unique(attractor_list);
        basin_matrix_list{x(1),x(2)} = attractor_list; % save to for efficiency
    end
end


function [region_list,changed] = merge_regions(perimeter_list,region_list,D,D1,D2)

changed = false;
for i=1:length(perimeter_list)
    I = find((D(i,:)<2.5) & (D1(i,:)<1.5) & (D2(i,:)<1.5)); %find all neighbours;
    if sum(region_list(I)>min(region_list(I)))>0
        region_list(I) = min(region_list(I));
        changed = true;
    end
end


end
