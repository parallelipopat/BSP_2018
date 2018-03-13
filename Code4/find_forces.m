function [forces, potential_energy] = find_forces(num_particles, epsilon, sigma, coordinates, length_cube, neighbours_list, num_neighbours_list)
    potential_energy = 0;
    sigma_6 = sigma^6;
    forces = zeros(3*num_particles, 1);
    
    for i = 1:3:3*num_particles
        i_index = (i+2)/3;
        for neighbour = 1:num_neighbours_list(i_index)
            j = neighbours_list(i_index, neighbour);
            if (j > i)
                diff_r = coordinates(i:i+2) - coordinates(j:j+2);
                diff_r = diff_r - length_cube*round(diff_r/length_cube);
                dist_r_2 = sum(diff_r.^2);
                dist_r_6 = dist_r_2^3; dist_r_8 = dist_r_6 * dist_r_2;
                dist_r_10 = dist_r_8 * dist_r_2;
                
                force_factor = (sigma_6/dist_r_8)*(sigma_6/dist_r_6 - 0.5);
                forces(i:i+2) = forces(i:i+2) + force_factor*diff_r;
                forces(j:j+2) = forces(j:j+2) - force_factor*diff_r;
                               
                potential_energy = potential_energy + (sigma_6/dist_r_6)*((sigma_6/dist_r_6) - 1);
            end
        end
    end
    forces = 48*epsilon*forces;
    potential_energy = 4*epsilon*potential_energy;
end