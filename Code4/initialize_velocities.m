function velocities = initialize_velocities(num_particles, mass, kB, temperature)
    seed = 7;
    rng(seed);
    
    velocities_3d = normrnd(0, (kB*temperature/mass)^0.5, num_particles, 3);
    velocities_3d = velocities_3d - sum(velocities_3d)/num_particles;
    velocities = [velocities_3d(:,1)'; velocities_3d(:,2')'; velocities_3d(:,3)'];
    velocities = velocities(:);
end
