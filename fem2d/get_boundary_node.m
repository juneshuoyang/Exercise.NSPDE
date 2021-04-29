%Function which tells us when a point is on the boundary of the domain.
%If not on boundary it returns 0. If the boundary is Dirichlet
%returns 1, if the boundary is Neumann returns 2.

function boundary_node = get_boundary_node(coords)

%Tolerance
eps = 1.0E-7;

 if (abs(coords(1)-0.0)<eps | abs(coords(2)-0.0) <eps | abs(coords(2)-1.0) <eps ...
         | abs(coords(1)-1.0) < eps)
     %Dirichlet Boundary
     boundary_node = 1;
 else
     %Not on boundary
     boundary_node = 0;
 end

