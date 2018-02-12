function [coordinates, length_cube] = initialize_cube(num_particles, density)
    %Rahman's parameters: 864, 1.374 g cm^-3
    
    mass_one_argon_atom = 39.95 * 1.6747 * (10^-24); % grams
    sigma_argon = 3.4 * (10^-8); % cm
    length_cube = ((num_particles * mass_one_argon_atom) / (density * sigma_argon^3))^(1/3); % in terms of sigma
    
    seed = 7;
    rng(seed);
    
    coordinates = length_cube.*rand(num_particles, 3);
    scatter3(coordinates(:, 1), coordinates(:, 2), coordinates(:, 3), 'filled');
end