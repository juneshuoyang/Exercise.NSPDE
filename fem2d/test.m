limits = [0 1; 0 1];

number_points(1) = input('Please enter number nodes in x-direction ');
number_points(2) = input('Please enter number nodes in y-direction ');
ele_type = input('Please enter element type, 3 = triangle, 4 quadrilateral ');

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
