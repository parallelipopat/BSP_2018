function velocities = initialize_velocities(num_particles, temperature)
    mass_one_atom = 48;
    kB = 1;
    
    seed = 7;
    rng(seed);
    
    velocities = normrnd(0, (kB*temperature/mass_one_atom)^0.5, num_particles, 3);
end
