function [forces, potential_energy] = find_forces(num_particles, epsilon, sigma, coordinates, length_cube, neighbours_list, num_neighbours_list)
    potential_energy = 0;
    sigma_6 = sigma^6;
    forces = zeros(num_particles, 3);
    for i = 1:num_particles-1
        for neighbour = 1:num_neighbours_list(i)
            j = neighbours_list(i, neighbour);
            if (j > i)
                diff_r = coordinates(i,:) - coordinates(j,:);
                diff_r = diff_r - length_cube*round(diff_r/length_cube);
                dist_r_2 = sum(diff_r.^2);
                dist_r_6 = dist_r_2^3; dist_r_8 = dist_r_6 * dist_r_2;
                dist_r_12 = dist_r_6^2; dist_r_14 = dist_r_12 * dist_r_2;
                
                force_factor = (1/dist_r_14) - 0.5*(1/dist_r_8);
                forces(i,:) = forces(i,:) + force_factor*diff_r;
                forces(j,:) = forces(j,:) - force_factor*diff_r;
                potential_energy = potential_energy + (1/dist_r_6)*((1/dist_r_6) - 1);
            end
        end
    end
    forces = forces*48*epsilon;
    potential_energy = potential_energy*4*epsilon;
end