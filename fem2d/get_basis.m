%Function to return the values of the basis functions and their gradients
%on the reference element at position local_coords

function [basis_funcs,grad_basis_funcs] = get_basis(local_coords,nvert,no_points)

if (nvert==3) %Trangular mesh
    %Calculate basis functions - rows are basis number, columns quad points
    basis_funcs(1,:) = 1 - local_coords(1, :) - local_coords(2, :); %%% COMPLETE HERE
    basis_funcs(2,:) = local_coords(1, :); %%% COMPLETE HERE
    basis_funcs(3,:) = local_coords(2, :); %%% COMPLETE HERE
    
    %Calculate grad of basis function - first column basis number, 2nd x/y
    %derivative, 3rd quad points
    
    %First basis
    grad_basis_funcs(1,1,1:no_points) = -1; %%% COMPLETE HERE
    grad_basis_funcs(1,2,1:no_points) = -1; %%% COMPLETE HERE
    %Second basis
    grad_basis_funcs(2,1,1:no_points) = 1; %%% COMPLETE HERE
    grad_basis_funcs(2,2,1:no_points) = 0; %%% COMPLETE HERE
    %Thirds basis
    grad_basis_funcs(3,1,1:no_points) = 0; %%% COMPLETE HERE
    grad_basis_funcs(3,2,1:no_points) = 1; %%% COMPLETE HERE
elseif (nvert==4) %Quadrilateral mesh
    %basis functions
    basis_funcs(1,:) = (local_coords(1,:).*local_coords(2,:)-local_coords(1,:)-...
        local_coords(2,:)+1)/4;
    basis_funcs(2,:) = (-local_coords(1,:).*local_coords(2,:)+local_coords(1,:)-...
        local_coords(2,:)+1)/4;
    basis_funcs(3,:) = (local_coords(1,:).*local_coords(2,:)+local_coords(1,:)+...
        local_coords(2,:)+1)/4;
    basis_funcs(4,:) = (-local_coords(1,:).*local_coords(2,:)-local_coords(1,:)+...
        local_coords(2,:)+1)/4;
    
    %grad of basis functions
    %First basis function
    grad_basis_funcs(1,1,1:no_points) = (local_coords(2,:)-1)/4;
    grad_basis_funcs(1,2,1:no_points) = (local_coords(1,:)-1)/4;
    %Second basis function
    grad_basis_funcs(2,1,1:no_points) = (-local_coords(2,:)+1)/4;
    grad_basis_funcs(2,2,1:no_points) = (-local_coords(1,:)-1)/4;
    %Third basis function
    grad_basis_funcs(3,1,1:no_points) = (local_coords(2,:)+1)/4;
    grad_basis_funcs(3,2,1:no_points) = (local_coords(1,:)+1)/4;
    %fourth basis function
    grad_basis_funcs(4,1,1:no_points) = (-local_coords(2,:)-1)/4;
    grad_basis_funcs(4,2,1:no_points) = (-local_coords(1,:)+1)/4;
end
end