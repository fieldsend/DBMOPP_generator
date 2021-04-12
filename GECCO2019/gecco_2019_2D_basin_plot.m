function [perimeter_list, optima_list, region_list, mode_matrix, basin_matrix, B] ...
    = gecco_2019_2D_basin_plot(dpp, n, function_flag, num_objs, min_v, max_v)
 
% function [perimeter_list, optima_list, region_list, mode_matrix, basin_matrix, B] ...
%    = gecco_2019_2D_basin_plot(dpp, n, function_flag, num_objs, min_v, max_v)   
%
% INPUTS
%
% dpp = distance-based point structure from the generator
% n = number of samples per axis (grid resolution)
% function_flag = (OPTIONAL) if argument value is 0, dpp is treated as a 
%     distance-based point structure instance from the generator. If the 
%     argument is 1, dpp is treated as a function name in a string, and 
%     will be invoked with feval assumingthe form y = f(x,num_objs). If the
%     argument value is 2, dpp is assumed to be a n by n by num_objs matrix 
%     holding precomputed values. Default argument value if not supplied is
%     0.
% num_objs = (OPTIONAL) used if function_flag is true. Number of objectives
%     in function argument. Will take from distance-based point structure
%     if not supplied.
% min_v = (OPTIONAL) used if function_flag is true. Minimum design space 
%     values (box constraint lower bound). Will use -1 if not supplied, as 
%     distance-based point structure assumed
% max_v = (OPTIONAL) used if function_flag is true. Maximum design space 
%     values (box constraint upper bound). Will use 1 if not supplied, as 
%     distance-based point structure assumed
%
% OUTPUTS
%
% perimeter_list = m by 2 matrix, each row holds indices corresponding to a
%                  cell in the n by n matrix of design locations, which is
%                  on the boundary (adjacent) of the dominance neutral
%                  region -- i.e., it is dominance neutral (ie. it has no
%                  dominanting neighbours), but it has at least one
%                  neighbour it dominates. It will therefore be a subset of
%                  the optima_list (including the possibility of being 
%                  equal).
% optima_list =  k by 2 matrix, each row holds a indices corresponding to a
%                  cell in the n by n matrix of design locations, which is
%                  a dominance neutral cell
% region_list = m by 1 vector. Holds a value representing the particular
%                  contiguous dominance neutral region the corresponding 
%                  element of perimeter_list lies on the perimeter of 
% mode_matrix = n by n matrix. A cell has a value of 0 if it is dominated
%                  by a neighbour. A value of 1 if it is mutually
%                  non-dominating with all neighbours, and avalue of 2
%                  otherwise (it is not dominated, but has at least one
%                  neighbour it dominates).
% basin_matrix = n by n matrix. Holds the number of different neutral 
%                  regions reachable from dominating paths downhill. A 0 
%                  value is therefore associated with a dominance neutral 
%                  cell.
% B = n by n matrix. dominance landscape matrix. A cell value of -1 
%                  indicates the cell is a dominance neutral optima, all 
%                  neighbours are mututally non-dominating, or dominated. 
%                  A cell value of 1 denotes
%                  that the cell is in the basin of a single dominance
%                  neutral region -- all dominating 'downhill' paths from
%                  the cell lead to the same contigous domiance neutral
%                  region. A value of 1 indicates the cell is on a basin
%                  boundary/saddle point, with multiple different neutral 
%                  regions reachable from dominating paths downhill.
%
% Generates various 2D basin plots for the given problem structue at the 
% grid resolution provided, from either a distance-based point structure
% instance from the generator, or another function (assumed to take a two 
% element design vector, and the number of objectives as an argument)
%
% Uses dits2 from Ian Nabney's Netlab toolbox
%
% Jonathan Fieldsend, University of Exeter, 2018, 2019, 2021
% See license information in package, available at 
% https://github.com/fieldsend/DBMOPP_generator

if exist('function_flag','var') == false
    function_flag = 0;
end
if (function_flag == 0)    
    x = linspace(-1,1,n);
    min_v = -1;
    max_v = 1;
    num_objs=dpp.num_objectives;
else
    x = linspace(min_v,max_v,n);
end

R = zeros(n,n,num_objs);

if (function_flag==2)
    R = dpp;
else
    % get qualities for each
    for i=1:n
        for j=1:n
            if (function_flag==1)
                t = feval(dpp,[x(i), x(j)],num_objs);
            else
                t = distance_points_problem([x(i), x(j)],dpp);
            end
            R(i,j,:) = t';
        end
    end
end

% now identify local optima regions
[perimeter_list, optima_list,region_list,mode_matrix,basin_matrix] = get_perimeter_and_optima(R,n);

figure; plot((perimeter_list(:,1)/(n))*(max_v-min_v)+min_v,(perimeter_list(:,2)/(n))*(max_v-min_v)+min_v,'k.');
axis([min_v max_v min_v max_v;]);
axis square;
figure; plot((optima_list(:,1)/(n))*(max_v-min_v)+min_v,(optima_list(:,2)/(n))*(max_v-min_v)+min_v,'r.');
axis([min_v max_v min_v max_v;]);
axis square;

figure; hold on;  imagesc(mode_matrix')
set(gca,'YDir','normal')

figure;hold on;  imagesc(basin_matrix')
set(gca,'YDir','normal')

B = basin_matrix-1;
B(B>1) = 1; 

figure; pcolor(x,x,B');%imagesc(B')
shading flat;
set(gca,'YDir','normal');
ylabel('$x_2$','Interpreter','latex');
xlabel('$x_1$','Interpreter','latex');
set(gca,'FontSize',20);
caxis([-1 1]); % in case a single basin of attraction


end

function x = dominates(a,b)
x = false;
if (sum(a<b)>0) && (sum(a<=b)==length(a))
    x = true;
end
end


function [perimeter_list, optima_list, region_list, mode_matrix, basin_matrix] = get_perimeter_and_optima(R,n)

perimeter_list = zeros(0,2);
optima_list = zeros(0,2); % each row will hold a coordinate of a dominance neutral cell
mode_matrix = zeros(n,n); % will be 1 in a corresponding dominance neutral 
% cell of R, 2 if direcly adjacent to a dominance neutral cell, zero otherwise

basin_matrix = zeros(n,n);
perimeter_matrix = zeros(n,n);

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
% optima list now contains index pairs (coordinates) of dominance neutral 
% cells, and mode_matrix cells are 1 in these locations, 0 in all others 

% determine neighbours of each cell
Neighbours=cell(n,n);
dominating_neighbours=cell(n,n);
for i=1:n
    for j=1:n
        indices = [];
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

% Could reorder the preceeding two blocks for efficiency...

%now determine perimeter list
for k=1:length(optima_list)
    % process each dominance neutral cell in turn
    i = optima_list(k,1);
    j = optima_list(k,2);
    %perimeter_flag = true;
    if length(Neighbours{i,j})~=8 % if a mode is on the boundary of the box constraint domain, it is automatically a perimeter member
         perimeter_list(end+1,:) = [i,j];
         mode_matrix(i,j) = 2;
    else
        for m=1:length(Neighbours{i,j})
            % if a neighbour of a neutral cell is not also a neutral cell, flag
            % as a perimeter cell, bounding the neutral region
            if (mode_matrix(Neighbours{i,j}(m,1),Neighbours{i,j}(m,2))==0)
                perimeter_list(end+1,:) = [i,j];
                mode_matrix(i,j) = 2;
                break; % don't need to process (i,j) any more, as already flagged on perimeter
            end
        end
    end
end

fprintf('iterate\n');
% iterative determine region -- worth checking if quicker way!

region_list = (1:length(perimeter_list))';
changed = true;
D = dist2(perimeter_list,perimeter_list);
D1 = dist2(perimeter_list(:,1),perimeter_list(:,1));
D2 = dist2(perimeter_list(:,2),perimeter_list(:,2));
%D(D==0) = inf; % set diagonal to infinity
step = 0;
% now iteratively merge together the neutral cells into disconnected
% regions
while(changed==true)
    [region_list,changed] = merge_regions(perimeter_list,region_list,D,D1,D2);
    step=step+1;
    fprintf('iteration %d, number of regions %d\n',step, length(unique(region_list)));
end
% region_list holds an (arbitary) numeric for each perimeter cell, based on the region it
% bounds

% label perimeter
for k=1:length(perimeter_list)
    perimeter_matrix(perimeter_list(k,1),perimeter_list(k,2)) = region_list(k);
end

basin_matrix_list = cell(n,n);    
% now compute and store the number of unique basins each cell flows to
for i=1:n
    fprintf('matrix row %d processing... \n',i);
    for j=1:n
        if (mode_matrix(i,j)==0) % if not a mode, will flow somewhere as has at least one dominating neighbour
            [attractor_list, basin_matrix_list] = recursive_attractor([i j], dominating_neighbours, mode_matrix,perimeter_matrix,basin_matrix_list);
            a = unique(attractor_list);
            basin_matrix(i,j) = length(a);
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
% process each element of the perimeter list
for i=1:length(perimeter_list)
    I = find((D(i,:)<2.5) & (D1(i,:)<1.5) & (D2(i,:)<1.5)); %find all neighbours of ith in design space;
    if sum( region_list(I) > min(region_list(I)) )>0
        region_list(I) = min(region_list(I));
        changed = true;
    end
end


end
