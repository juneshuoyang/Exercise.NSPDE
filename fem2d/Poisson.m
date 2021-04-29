clear all
limits = [0 1; 0 1];

%User input
number_points(1) = input('Please enter number nodes in x-direction ');
number_points(2) = input('Please enter number nodes in y-direction ');
ele_type = input('Please enter element type, 3 = triangle, 4 quadrilateral ');

no_vert = ele_type;
%Number Gauss points in each direction
n = 4;

%Set the colourmap
cm = colormap(jet(256));

%Get the quadrature rule for the reference element - Gauss quadrature is used
[quad_points,quad_weights,no_points] = gauss_quadrature(n,ele_type);

%Call mesh generator
[no_eles,no_nodes,coords,ele_connect,boundary] = mesh_generator(limits,number_points,...
    ele_type);

%plot mesh
figure
hold on
for k=1:no_eles
    %Store in "ele_nodes" the node numbers for the element k
    ele_nodes = ele_connect(k, :);%%% COMPLETE HERE
    
    %Store in "local_coords" the coordinates of the element nodes
	%%% COMPLETE HERE
    local_coords = coords(:, ele_nodes);
    
    % Plot
	plot([local_coords(1,:) local_coords(1,1)],[local_coords(2,:) local_coords(2,1)])

end

%Initialize matrix A as sparse matrix and rhs
A = sparse(no_nodes, no_nodes);
rhs = zeros(no_nodes, 1);

%Loop over the elements
for k=1:no_eles
    %Store in "ele_nodes" the node numbers for the element k
    %%% COMPLETE HERE
    ele_nodes = ele_connect(k, :);
    
    %Store in "local_coords" the coordinates of the element nodes
    %%% COMPLETE HERE
    local_coords = coords(:, ele_nodes);

    %Caculate the element Jacobian
    %Here we assume an affine mapping from reference ele to global ele
    % x = Jac*x_ref+bvec
    [jacobi_mat,jacobian,bvec] = get_jacobian(local_coords,no_vert);
      
    
    %Get basis functions and gradients in local coordinates... REF
    [basis_funcs,grad_basis_funcs] = get_basis(quad_points,no_vert,no_points);
    
    %Map local quadrature points to global coordinates
    global_points = map_ref(quad_points,jacobi_mat,bvec,no_points);
    
    %Change gradients to global coordinate system
    inv_jacobi = inv(jacobi_mat');
    tmp=zeros(2,no_points);
    for i=1:no_vert
        tmp(:,:)=grad_basis_funcs(i,:,:);
        global_grad(i,:,:) = inv_jacobi*tmp;
    end
			
    
    %Setup the local stiffness matrix
    %Initialize to zero
    local_stiff = zeros(no_vert,no_vert);
    local_rhs = zeros(no_vert,1);
    %Loop over the quadrature points (could be optimised)
    for m=1:no_points
        %calculate source term at quadrature point
        source_term = source_data(global_points(:, m)); %%% COMPLETE HERE
        %Loop over test basis functions
        for i=1:no_vert
            %Loop over trial basis functions
            for j=1:no_vert
                local_stiff(i,j) = local_stiff(i,j) + global_grad(i,:,m) * global_grad(j,:,m)' * quad_weights(m) * jacobian; %%% COMPLETE HERE
            end
        local_rhs(i) = local_rhs(i) + source_term * basis_funcs(i, m) * quad_weights(m) *  jacobian; %%% COMPLETE HERE
        end
    end
    
    %Distribute local stiffness to global stiffness matrix
    %%% COMPLETE HERE
					
    %Distribute local rhs to global rhs
	%%% COMPLETE HERE
    
    A(ele_nodes, ele_nodes) = A(ele_nodes, ele_nodes) + local_stiff;
    rhs(ele_nodes) = rhs(ele_nodes) + local_rhs; 

end
			
%Loop over nodes and modify global matrix/RHS if a Dirichlet node
% inserting value using boundary_data function
%%% COMPLETE HERE		

b = 1 : no_nodes;
b = b .* boundary;
b = b(b~=0);

A(b, :) = 0; A(:, b) = 0; rhs(b) = 0;
A(b, b) = speye(length(b), length(b));

%Solve equation for unknowns U
U = A\rhs;

%Plot solution
 figure
 hold on
 %One element at a time
 for k=1:no_eles
     ele_nodes=ele_connect(k,:);
        coords1 = [coords(1,ele_nodes) coords(1,ele_nodes(1))];
        coords2 = [coords(2,ele_nodes) coords(2,ele_nodes(1))];
        z = [U(ele_nodes); U(ele_nodes(1))]';

      %Plot with a single colour - green in this case
      fill3(coords1,coords2,z,'g')
      
      %Plot with interpolated colours
      %fill3(coords1,coords2,z,z,'EdgeColor','none')

 end
 
 view([-1,-1,1]);
