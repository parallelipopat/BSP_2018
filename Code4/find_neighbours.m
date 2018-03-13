function [neighbours_list,num_neighbours_list] = find_neighbours(num_particles, coordinates, length_cube, r_cutoff)
    neighbours_list = zeros(num_particles, num_particles -1);
    num_neighbours_list = zeros(num_particles, 1);
    r_cutoff_2 = r_cutoff*r_cutoff;
    for i = 1:3:3*num_particles
        num_neighbours_i = 0; 
        i_index = (i+2)/3;
        for j = i+3:3:3*num_particles
            j_index = (j+2)/3;
            diff_r = coordinates(i:i+2) - coordinates(j:j+2);
            diff_r = diff_r - length_cube*round(diff_r/length_cube);
            diff_r_2 = sum(diff_r.^2);
            if diff_r_2 < r_cutoff_2
                num_neighbours_i = num_neighbours_i + 1;
                neighbours_list(i_index, num_neighbours_i) = j;
                num_neighbours_j = num_neighbours_list(j_index);
                num_neighbours_j = num_neighbours_j + 1;
                neighbours_list(j_index, num_neighbours_j) = i;
                num_neighbours_list(j_index) = num_neighbours_j;
            end
        end
        num_neighbours_list(i_index) = num_neighbours_i;
    end
    neighbours_list = neighbours_list(:,1:max(num_neighbours_list));
end