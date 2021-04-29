function [basis_funcs, grad_basis_funcs] = get_basis_1d(ref_points)
    
    basis_funcs(1, :) = - 1 / 2 * ref_points + 1 / 2;
    basis_funcs(2, :) = 1 / 2 * ref_points + 1 / 2;
    
    n = size(ref_points, 2);
    
    %First basis
    grad_basis_funcs(1, 1 : n) = - 1 / 2; 
    
    %Second basis
    grad_basis_funcs(2, 1 : n) = 1 / 2; 
    
end
    