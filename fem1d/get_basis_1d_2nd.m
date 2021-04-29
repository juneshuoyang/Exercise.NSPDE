function [basis_funcs, grad_basis_funcs] = get_basis_1d_2nd(ref_points)
    
    basis_funcs(1, :) = 1 / 2 * (ref_points - 1) .* ref_points;
    basis_funcs(2, :) = - (ref_points + 1) .* (ref_points - 1);
    basis_funcs(3, :) = 1 / 2 * (ref_points + 1) .* ref_points;
    
    %First basis
    grad_basis_funcs(1, :) = ref_points - 1 / 2; 
    
    %Second basis
    grad_basis_funcs(2, :) = - 2 * ref_points; 
    
    %Third basis
    grad_basis_funcs(3, :) = ref_points + 1 / 2; 
    
end
    