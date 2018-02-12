function forces = LJ_Force(coordinates, length_cube)
    forces = zeros(size(coordinates));
    num_particles = size(coordinates, 1);
    
    for i = 1:num_particles-1
        for j = i+1:num_particles
            delta_r = coordinates(i,:) - coordinates(j,:);
            adjusted_delta_r = periodic_boundary_correction(delta_r, length_cube);
            inverse_delta_r = 1/norm(adjusted_delta_r);
            force_ij = (inverse_delta_r^8 * (inverse_delta_r^6 - 0.5)) * adjusted_delta_r;
            
            forces(i,:) = forces(i,:) + force_ij;
            forces(j,:) = forces(j,:) - force_ij;
        end
    end
    
    forces = forces * 48;
end