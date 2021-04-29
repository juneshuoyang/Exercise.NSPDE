function [jacobi_mat, jacobian, bvec] = get_jacobian_1d(ele_coords)
 
    jacobi_mat = (ele_coords(2) - ele_coords(1)) / 2;
    jacobian = det(jacobi_mat);
    bvec = (ele_coords(2) + ele_coords(1)) / 2;
    
end