function velocities = initialize_velocities(num_particles, temperature)
    [~, ~, ~, mass, ~, kB, ~, ~, ~] = initialise_params();
    
    seed = 7;
    rng(seed);
    
    velocities = normrnd(0, (kB*temperature/mass)^0.5, num_particles, 3);
    velocities = velocities - sum(velocities)/num_particles;
end
