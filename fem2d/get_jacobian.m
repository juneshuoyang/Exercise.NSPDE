%Function to return the Jacobi matrix and jacobian of the element mapping
%For quads we assume a parallelogram

function [jacobi_mat,jacobian,bvec] = get_jacobian(ele_coords,nvert)

if (nvert==3) %Triangular Mesh - Reference element is [(0,1),(1,0),(0,1)]
	%%% COMPLETE HERE
    jacobi_mat = [(ele_coords(1,2)-ele_coords(1,1)) (ele_coords(1,3)-ele_coords(1,1));...
    (ele_coords(2,2)-ele_coords(2,1)) (ele_coords(2,3)-ele_coords(2,1))];

    jacobian = det(jacobi_mat);
    bvec = [ele_coords(1,1); ele_coords(2,1)];
		
elseif (nvert==4) %Quadrialteral Mesh - Reference element is [(-1,-1),(1,-1),(1,1),(-1,1)]
    %Here we assume the quadrilateral is a parallelogram and hence the
    %mapping is affine
    jacobi_mat = [(ele_coords(1,2)-ele_coords(1,1))/2 (ele_coords(1,4)-ele_coords(1,1))/2;...
        (ele_coords(2,2)-ele_coords(2,1))/2 (ele_coords(2,4)-ele_coords(2,1))/2];
    
    jacobian = det(jacobi_mat);
    bvec = [(ele_coords(1,2)+ele_coords(1,4))/2;(ele_coords(2,2)+ele_coords(2,4))/2];
end
end
