function velocities = initialize_velocities(num_particles, mass, kB, temperature)
    % Seeding to maintian reproducibility
    seed = 7;
    rng(seed);
    
    % Generate velocities using Maxwell distribution, and remove linear momentum
    velocities_3d = normrnd(0, (kB*temperature/mass)^0.5, num_particles, 3);
    velocities_3d = velocities_3d - sum(velocities_3d)/num_particles;
    %Arrange as a column vector
    velocities = [velocities_3d(:,1)'; velocities_3d(:,2')'; velocities_3d(:,3)'];
    velocities = velocities(:);
end
