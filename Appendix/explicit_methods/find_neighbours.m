function [neighbours_list,num_neighbours_list] = ...
          find_neighbours(num_particles, coordinates, length_cube, r_cutoff)
    neighbours_list = zeros(num_particles, num_particles -1);
    num_neighbours_list = zeros(num_particles, 1);
    r_cutoff_2 = r_cutoff*r_cutoff;
    for i = 1:num_particles-1
        num_neighbours_i = 0; 
        for j = i+1:num_particles
            diff_r = coordinates(i,:) - coordinates(j, :);
            diff_r = diff_r - length_cube*round(diff_r/length_cube);
            dist_r_2 = sum(diff_r.^2);
            if dist_r_2 < r_cutoff_2
                num_neighbours_i = num_neighbours_i + 1;
                neighbours_list(i, num_neighbours_i) = j;
                num_neighbours_j = num_neighbours_list(j);
                num_neighbours_j = num_neighbours_j + 1;
                neighbours_list(j, num_neighbours_j) = i;
                num_neighbours_list(j) = num_neighbours_j;
            end
        end
        num_neighbours_list(i) = num_neighbours_i;
    end
    % Re-adjust size of neighbours_list to save memory
    neighbours_list = neighbours_list(:,1:max(num_neighbours_list));
end