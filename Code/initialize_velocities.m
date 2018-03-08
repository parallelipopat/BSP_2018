function velocities = initialize_velocities(num_particles, temperature)
    [num_particles, epsilon, sigma, mass, density, kB, temperature, h] = initialise_params()
    
    seed = 7;
    rng(seed);
    
    velocities = normrnd(0, (kB*temperature/mass)^0.5, num_particles, 3);
end
