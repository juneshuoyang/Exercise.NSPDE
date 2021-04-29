% Exercise 2.1

a = 0;
b = pi / 2; % 1;
limits = [a b];

nGaus = 4;
[quad_points, quad_weights] = gauss_1d(nGaus);
[basis_funcs, grad_basis_funcs] = get_basis_1d(quad_points);

nG_ea = 6; % error analysis
[q_ea, w_ea] = gauss_1d(nG_ea);
[basis_funcs_ea, grad_basis_funcs_ea] = get_basis_1d(q_ea);

fprintf('nELm           L2             H1\n');

for l = 1 : 5

    nElm = 10 * 2 ^ (l - 1);
    nNod = nElm + 1;
    coords = linspace(a, b, nNod);

    A = zeros(nNod, nNod);
    rhs = zeros(nNod, 1);

    for k = 1 : nElm

        ele_nodes = [k k + 1];
        local_coords = coords(ele_nodes);

        [jacobi_mat, jacobian, bvec] = get_jacobian_1d(local_coords);   
        global_points = map_ref_1d(quad_points, jacobi_mat, bvec);

        inv_jacobi = 1 / jacobi_mat; % jaboci_mat in 1D is a scalar
        global_grad = inv_jacobi * grad_basis_funcs;     

        local_stiff = zeros(2, 2);
        local_rhs = zeros(2, 1);

        for m = 1 : nGaus
            source_term = source_data_1d(global_points(m));
            for i = 1 : 2
                for j = 1 : 2
                    local_stiff(i, j) = local_stiff(i, j) + global_grad(i, m) * global_grad(j, m) * quad_weights(m) * jacobian;
                end
            local_rhs(i) = local_rhs(i) + source_term * basis_funcs(i, m) * quad_weights(m) *  jacobian; 
            end
        end

        A(ele_nodes, ele_nodes) = A(ele_nodes, ele_nodes) + local_stiff;
        rhs(ele_nodes) = rhs(ele_nodes) + local_rhs; 

    end 

    rhs(1) = exact_data_1d(a); rhs(2) = rhs(2) - exact_data_1d(a) * A(2, 1);
    A(1, :) = 0; A(:, 1) = 0; 
    A(1, 1) = 1;

    rhs(end) = exact_data_1d(b); rhs(end - 1) = rhs(end - 1) - exact_data_1d(b) * A(end - 1, end);
    A(end, :) = 0; A(:, end) = 0; 
    A(end, end) = 1;

    U = A \ rhs;

    plot(coords, U, '.'); 
    title('Numerical solution of $-u^{\prime\prime} = \sin(x)$ on [0, $\pi/2$]', 'Interpreter', 'latex')

    L2 = 0;
    H1 = 0;
    for k = 1 : nElm

        ele_nodes = [k k + 1];
        local_coords = coords(ele_nodes);

        [jacobi_mat, jacobian, bvec] = get_jacobian_1d(local_coords);   
        global_points = map_ref_1d(q_ea, jacobi_mat, bvec);
        
        inv_jacobi = 1 / jacobi_mat;
        global_grad_ea = inv_jacobi * grad_basis_funcs_ea;

        f_err_elem = exact_data_1d(global_points) - U(k) * basis_funcs_ea(1, :) - U(k + 1) * basis_funcs_ea(2, :);
        g_err_elem = prime_data_1d(global_points) - U(k) * global_grad_ea(1, :) - U(k + 1) * global_grad_ea(2, :);
        L2 = L2 + w_ea * (f_err_elem .* f_err_elem)' * jacobian;
        H1 = H1 + w_ea * (f_err_elem .* f_err_elem)' * jacobian + w_ea * (g_err_elem .* g_err_elem)' * jacobian;
        
    end

    L2 = sqrt(L2); H1 = sqrt(H1);
    Iu = exact_data_1d(coords)';
    %Enrg = (Iu - U)' * A * (Iu - U); 
    %disp(Iu - U);
    fprintf('%d   %10.10f   %10.10f\n', nElm, L2, H1);
    
end