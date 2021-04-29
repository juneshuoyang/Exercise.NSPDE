%Simple rectangular mesh generator - Triangular or quadrilateral meshes
%Returns no_eles, no_nodes, the coords of each node, the ele connectivity
%of each element and whether each node is a boundary node or not, plus also
%return no_boundary_edges and edge_info - nodes of edges and edge type -
%1 = Dirichlet, 2 = Neumann 

function [no_eles,no_nodes,coords,ele_connect,boundary,no_boundary_edges,...
    edges] = mesh_generator(limits,number_points,ele_type)

%Setup node coordinates
%Ordered from left to right, then down to up
count = 1;
h1 = (limits(1,2)-limits(1,1))/(number_points(1)-1);
h2 = (limits(2,2)-limits(2,1))/(number_points(2)-1);
for j=1:number_points(2)
    for i=1:number_points(1)
        coords(1,count) = limits(1,1)+(i-1)*h1;
        coords(2,count) = limits(2,1)+(j-1)*h2;
        %Check if on the boundary or not
        boundary(count) = get_boundary_node(coords(:,count));
        
        count = count+1;
    end
end

no_nodes = count-1;

%Setup element connectivity and boundary edges
bottom_left = 1;
bottom_right = 2;
top_left = 1+number_points(1);
top_right = 2+number_points(1);
count = 1;
count_edges = 1;
for j=1:number_points(2)-1
    for i=1:number_points(1)-1
        if (ele_type==3) %Triangular mesh
            %First element
            ele_connect(count,1) = bottom_left;
            ele_connect(count,2) = bottom_right;
            ele_connect(count,3) = top_left;
            
            %Setup edges by looping around element edges and testing
            %whether on boundary or not
            if (get_boundary_node((coords(:,bottom_left)+coords(:,bottom_right))/2)~=0)
                %First entry is boundary type, 2nd and 3rd entries are nodes,
                %4th is local edge number,5th is element which edge belongs
                %to
                edges(count_edges,1) = ...
                    get_boundary_node((coords(:,bottom_left)+coords(:,bottom_right))/2);
                edges(count_edges,2:4) = [bottom_left, bottom_right, 1];
                edges(count_edges,5) = count;
                count_edges = count_edges+1;
            end
            if (get_boundary_node((coords(:,bottom_right)+coords(:,top_left))/2)~=0)
                %First entry is boundary type, 2nd and 3rd entries are nodes,
                %4th is local edge number,5th is element which edge belongs
                %to
                edges(count_edges,1) = ...
                    get_boundary_node((coords(:,bottom_right)+coords(:,top_left))/2);
                edges(count_edges,2:4) = [bottom_right, top_left, 2];
                edges(count_edges,5) = count;
                count_edges = count_edges+1;
            end
            if (get_boundary_node((coords(:,top_left)+coords(:,bottom_left))/2)~=0)
                %First entry is boundary type, 2nd and 3rd entries are nodes,
                %4th is local edge number,5th is element which edge belongs
                %to
                edges(count_edges,1) = ...
                    get_boundary_node((coords(:,top_left)+coords(:,bottom_left))/2);
                edges(count_edges,2:4) = [top_left, bottom_left, 3];
                edges(count_edges,5) = count;
                count_edges = count_edges+1;
            end
            
            count = count+1;
            %Second element
            ele_connect(count,1) = top_right;
            ele_connect(count,2) = top_left;
            ele_connect(count,3) = bottom_right;
            
            %Setup boundary edges
            if (get_boundary_node((coords(:,top_right)+coords(:,top_left))/2)~=0)
                %First entry is boundary type, 2nd and 3rd entries are nodes,
                %4th is local edge number,5th is element which edge belongs
                %to
                edges(count_edges,1) = ...
                    get_boundary_node((coords(:,top_right)+coords(:,top_left))/2);
                edges(count_edges,2:4) = [top_right, top_left, 1];
                edges(count_edges,5) = count;
                count_edges = count_edges+1;
            end
            if (get_boundary_node((coords(:,top_left)+coords(:,bottom_right))/2)~=0)
                %First entry is boundary type, 2nd and 3rd entries are nodes,
                %4th is local edge number,5th is element which edge belongs
                %to
                edges(count_edges,1) = ...
                    get_boundary_node((coords(:,top_left)+coords(:,bottom_right))/2);
                edges(count_edges,2:4) = [top_left, bottom_right, 2];
                edges(count_edges,5) = count;
                count_edges = count_edges+1;
            end
            if (get_boundary_node((coords(:,bottom_right)+coords(:,top_right))/2)~=0)
                %First entry is boundary type, 2nd and 3rd entries are nodes,
                %4th is local edge number,5th is element which edge belongs
                %to
                edges(count_edges,1) = ...
                    get_boundary_node((coords(:,bottom_right)+coords(:,top_right))/2);
                edges(count_edges,2:4) = [bottom_right, top_right, 3];
                edges(count_edges,5) = count;
                count_edges = count_edges+1;
            end
            
            count = count+1;
        elseif (ele_type==4) %Quadrilateral mesh
            ele_connect(count,1) = bottom_left;
            ele_connect(count,2) = bottom_right;
            ele_connect(count,3) = top_right;
            ele_connect(count,4) = top_left;
            
            %Setup boundary edges
            if (get_boundary_node((coords(:,bottom_left)+coords(:,bottom_right))/2)~=0)
                %First entry is boundary type, 2nd and 3rd entries are nodes,
                %4th is local edge number,5th is element which edge belongs
                %to
                edges(count_edges,1) = ...
                    get_boundary_node((coords(:,bottom_left)+coords(:,bottom_right))/2);
                edges(count_edges,2:4) = [bottom_left, bottom_right, 1];
                edges(count_edges,5) = count;
                count_edges = count_edges+1;
            end
            if (get_boundary_node((coords(:,bottom_right)+coords(:,top_right))/2)~=0)
                %First entry is boundary type, 2nd and 3rd entries are nodes,
                %4th is local edge number,5th is element which edge belongs
                %to
                edges(count_edges,1) = ...
                    get_boundary_node((coords(:,bottom_right)+coords(:,top_right))/2);
                edges(count_edges,2:4) = [bottom_right, top_right, 2];
                edges(count_edges,5) = count;
                count_edges = count_edges+1;
            end
            if (get_boundary_node((coords(:,top_right)+coords(:,top_left))/2)~=0)
                %First entry is boundary type, 2nd and 3rd entries are nodes,
                %4th is local edge number,5th is element which edge belongs
                %to
                edges(count_edges,1) = ...
                    get_boundary_node((coords(:,top_right)+coords(:,top_left))/2);
                edges(count_edges,2:4) = [top_right, top_left, 3];
                edges(count_edges,5) = count;
                count_edges = count_edges+1;
            end
            if (get_boundary_node((coords(:,top_left)+coords(:,bottom_left))/2)~=0)
                %First entry is boundary type, 2nd and 3rd entries are nodes,
                %4th is local edge number,5th is element which edge belongs
                %to
                edges(count_edges,1) = ...
                    get_boundary_node((coords(:,top_left)+coords(:,bottom_left))/2);
                edges(count_edges,2:4) = [top_left, bottom_left, 4];
                edges(count_edges,5) = count;
                count_edges = count_edges+1;
            end
            
            count = count+1;
        end
        
        %Update nodes
        bottom_right = bottom_right+1;
        bottom_left = bottom_left+1;
        top_left = top_left+1;
        top_right = top_right+1;
    end
    %Go to beginning of next row
    bottom_right = bottom_right+1;
    bottom_left = bottom_left+1;
    top_left = top_left+1;
    top_right = top_right+1;
end

no_eles = count-1;
no_boundary_edges = count_edges-1;

end