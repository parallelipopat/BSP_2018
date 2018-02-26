function [coordinates, length_cube] = initialize_cube(num_particles, density)
    length_cube = (num_particles / density)^(1/3);
    
    seed = 7;
    rng(seed);
    
    coordinates = length_cube.*rand(num_particles, 3);
    scatter3(coordinates(:, 1), coordinates(:, 2), coordinates(:, 3), 'filled');
end