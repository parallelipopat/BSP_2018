function forces = LJ_Force(coordinates, length_cube)
    [num_particles, epsilon, sigma, ~, ~, ~, ~, ~, ~] = initialise_params();
    dimension = size(coordinates, 2);
    forces = zeros(num_particles, dimension);
    
    %d_forces = zeros(dimension, dimension, num_particles);

    for i = 1:num_particles-1
        for j = i+1:num_particles
            delta_r = coordinates(i,:) - coordinates(j,:);
            adjusted_delta_r = periodic_boundary_correction(delta_r, length_cube);
            
            inverse_delta_r = 1/norm(adjusted_delta_r);
            %force_ij = (inverse_delta_r^8 * (inverse_delta_r^6 - 0.5)) * adjusted_delta_r;
            force_ij = epsilon*((sigma*inverse_delta_r)^13 - 0.5*(sigma*inverse_delta_r)^7) * (inverse_delta_r) * adjusted_delta_r;
            forces(i,:) = forces(i,:) + force_ij;
            forces(j,:) = forces(j,:) - force_ij;

            %d_force_common_term = -14*inverse_delta_r^16 + 4*inverse_delta_r^10;
            %d_force_ij = d_force_common_term*(adjusted_delta_r'*adjusted_delta_r) + force_ij*eye(dimension);
            %d_forces(:,:,i) = d_forces(:,:,i) + d_force_ij;
            %d_forces(:,:,j) = d_forces(:,:,j) - d_force_ij;
        end
    end
    
    forces = forces * 48;
    %d_forces = d_forces * 48;
end