function velocities = initialize_velocities(num_particles, mass, kB, temperature)
    seed = 7;
    rng(seed);
    
    velocities = normrnd(0, (kB*temperature/mass)^0.5, num_particles, 3);
    velocities = velocities - sum(velocities)/num_particles;
end
