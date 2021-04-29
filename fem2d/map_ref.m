%Function which maps points in reference ele to global ele
%Assumes affine mapping x = jacobi_mat*x_ref+bvec

function global_points = map_ref(local_points,jacobi_mat,bvec,no_points)

%Multiply by jacobi_mat
%%% COMPLETE HERE

global_points = jacobi_mat * local_points;


%Add bvec
%%% COMPLETE HERE

global_points(:, 1 : no_points) = global_points(:, 1 : no_points) + bvec;

end
