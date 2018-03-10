function forces = LJ_Force_new(coordinates, length_cube)
    [num_particles, epsilon, sigma, ~, ~, ~, ~, ~, ~] = initialise_params();
    dimension = size(coordinates, 2);
    forces = zeros(num_particles, dimension);
    
    dim = 1;
    dim_coord = coordinates(:,dim);
    repeated_dim = repmat(dim_coord,1,num_particles);
    diff_dim = tril(transpose(repeated_dim) -  repeated_dim, -1);
    adjusted_diff_dim_x = arrayfun(@(z) periodic_boundary_correction_1d(z, length_cube), diff_dim);
    
    dim = 2;
    dim_coord = coordinates(:,dim);
    repeated_dim = repmat(dim_coord,1,num_particles);
    diff_dim = tril(transpose(repeated_dim) -  repeated_dim, -1);
    adjusted_diff_dim_y = arrayfun(@(z) periodic_boundary_correction_1d(z, length_cube), diff_dim);

    dim = 3;
    dim_coord = coordinates(:,dim);
    repeated_dim = repmat(dim_coord,1,num_particles);
    diff_dim = tril(transpose(repeated_dim) -  repeated_dim, -1);
    adjusted_diff_dim_z = arrayfun(@(z) periodic_boundary_correction_1d(z, length_cube), diff_dim);
    
    diff_sum = diff_sum + adjusted_diff_dim_x.*adjusted_diff_dim_x + adjusted_diff_dim_y.*adjusted_diff_dim_y  + adjusted_diff_dim_z.*adjusted_diff_dim_z;
    
    inv_diff_sum = 1./diff_sum;
    
    for i = 1:num_particles-1
        for j = i+1:num_particles
            squared_inverse_delta_r = inv_diff_sum(j, i);
            force_factor = (squared_inverse_delta_r^4 * (squared_inverse_delta_r^3 - 0.5)); 
            force_ij(1) = force_factor * adjusted_diff_dim_x(i, j);
            force_ij(2) = force_factor * adjusted_diff_dim_y(i, j);
            force_ij(3) = force_factor * adjusted_diff_dim_z(i, j);
            
            forces(i,:) = forces(i,:) + force_ij;
            forces(j,:) = forces(j,:) - force_ij;
        end
    end
    
    forces = forces * 48;
end
