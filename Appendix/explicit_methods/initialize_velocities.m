function velocities = initialize_velocities(num_particles, mass, kB, temperature)
    % Seeding for reproducibility
    seed = 7;
    rng(seed);
    
    % Generate velocities using Maxwell distribution, and remove linear momentum
    velocities = normrnd(0, (kB*temperature/mass)^0.5, num_particles, 3);
    velocities = velocities - sum(velocities)/num_particles;
end
