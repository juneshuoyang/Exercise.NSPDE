%Function to get 2D Gauss quadrature points and weights for reference
%triangle or reference quadrilateral

function [quad_points,quad_weights,no_points] = gauss_quadrature(poly_deg,ele_type)

%Calculate 1D Gaussian quadrature on interval [-1,1]
[points_1d, weights_1d] = gauss_1d(poly_deg);


%Tensor product the gauss points to the reference square
count = 1;

for j=1:poly_deg
    for i=1:poly_deg
        quad_points(1,count) = points_1d(i);
        quad_points(2,count) = points_1d(j);
        quad_weights(count) = weights_1d(i)*weights_1d(j);
        count = count+1;
    end
end

%Ref Triangle -s (-1,-1)->(1,-1)->(-1,1)
if (ele_type==3)  %Triangle - map points and weights from square to triangle
    quad_weights = (1-quad_points(2,:)).*quad_weights/2;
    quad_points(1,:) = (1-quad_points(2,:)).*(quad_points(1,:)+1)/2-1;
    
        %Now map to reference triangle (0,0)->(1,0)->(0,1)
    quad_weights = 0.25*quad_weights;
    quad_points = (quad_points+1)/2;
    	% patch if poly_deg = 1
	if (poly_deg == 1)
		quad_points = [1/3 ; 1/3];
    end
    
end

no_points = poly_deg^2;

end
        
        