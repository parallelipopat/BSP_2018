function coordinates = initialize_cube(num_particles, length_cube) 
    seed = 7;
    rng(seed);
    
    %a = -length_cube/2; b = -a;
    a = 0; b = length_cube;
    coordinates = (b-a).*rand(num_particles, 3)+ a;
    %scatter3(coordinates(:, 1), coordinates(:, 2), coordinates(:, 3), 'filled');
end