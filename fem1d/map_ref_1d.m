function global_points = map_ref_1d(ref_points, jacobi_mat, bvec)

    global_points = jacobi_mat * ref_points;

    n = size(ref_points, 2);
    global_points(1 : n) = global_points(1 : n) + bvec;

end
