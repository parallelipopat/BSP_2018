function [forces, potential_energy, jacobian_matrix] = find_forces(num_particles, epsilon, sigma, coordinates, length_cube, neighbours_list, num_neighbours_list)
    potential_energy = 0;
    sigma_6 = sigma^6;
    forces = zeros(3*num_particles, 1);
    jacobian_matrix = zeros(3*num_particles);
    
    for i = 1:3:num_particles
        for neighbour = 1:num_neighbours_list(i)
            j = neighbours_list(i, neighbour);
            if (j > i)
                diff_r = coordinates(i:i+2) - coordinates(j:j+2);
                diff_r = diff_r - length_cube*round(diff_r/length_cube);
                dist_r_2 = sum(diff_r.^2);
                dist_r_6 = dist_r_2^3; dist_r_8 = dist_r_6 * dist_r_2;
                dist_r_10 = dist_r_8 * dist_r_2;
                
                force_factor = (sigma_6/dist_r_8)*(sigma_6/dist_r_6 - 0.5);
                forces(i:i+2) = forces(i:i+2) + force_factor*diff_r;
                forces(j:j+2) = forces(j:j+2) - force_factor*diff_r;
                
                jacobian_factor = (sigma_6/dist_r_10)*((14*sigma_6/dist_r_6) - 4);
                jacobian_force_block = force_factor*eye(3) + jacobian_factor*(diff_r*diff_r');
                jacobian_matrix(i:i+2,i:i+2) = jacobian_matrix(i:i+2,i:i+2) + jacobian_force_block;
                jacobian_matrix(j:j+2,j:j+2) = jacobian_matrix(j:j+2,j:j+2) + jacobian_force_block;
                jacobian_matrix(i:i+2,j:j+2) = jacobian_matrix(i:i+2,j:j+2) - jacobian_force_block;
                jacobian_matrix(j:j+2,i:i+2) = jacobian_matrix(j:j+2,i:i+2) - jacobian_force_block;
                
                potential_energy = potential_energy + (sigma_6/dist_r_6)*((sigma_6/dist_r_6) - 1);
            end
        end
    end
    forces = 48*epsilon*forces;
    jacobian_matrix = 48*epsilon*jacobian_matrix;
    potential_energy = 4*epsilon*potential_energy;
end