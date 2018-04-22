clear variables;
tic;
[num_particles, epsilon, sigma, r_cutoff, mass, density, kB, temperature, h, N_e, N_f, N_s, N_n, beta, gamma] = initialize_params();
length_cube = find_cube_length(num_particles, density);
diff = length_cube/((2*num_particles)^(1/3));
coordinates = initialize_cube(num_particles, length_cube-diff);
velocities = initialize_velocities(num_particles, mass, kB, temperature);

[neighbours_list, num_neighbours_list] = find_neighbours(num_particles, coordinates, length_cube, r_cutoff);
[forces, ~] = find_forces(num_particles, epsilon, sigma, coordinates, length_cube, neighbours_list, num_neighbours_list);

energy = zeros(N_f/N_s, 3);
r_cutoff_6 = r_cutoff^6;
sigma_6 = sigma^6;
shifting_potential_term = 4*(sigma_6/r_cutoff_6)*((sigma_6/r_cutoff_6) - 1);

ax1 = subplot(1,2,1);
ax2 = subplot(1,2,2);
set(gcf, 'Position', get(0, 'Screensize'));

for i = 1:(N_e+N_f)
    if (mod(i, 25)==0)
        fprintf('%d\n',i);
    end
    
    acceleration = forces/mass;
    
    coordinates = coordinates + h*velocities + 0.5*h^2*acceleration;

    velocities = velocities + h*(1-gamma)*acceleration;
    [forces, potential_energy] = find_forces(num_particles, epsilon, sigma, coordinates, length_cube, neighbours_list, num_neighbours_list);
    acceleration = forces/mass;
    velocities = velocities + h*gamma*acceleration;
  
    if (i <= N_e)
       scaling_factor = (3*kB*num_particles*temperature)/(mass*sum(sum(velocities.^2,2)));
       velocities = velocities * sqrt(scaling_factor) ; 
    elseif (mod(i, N_s) == 0)
        energy_index = (i-N_e)/N_s;
        energy(energy_index, 1) = potential_energy/num_particles;
        energy(energy_index, 2) = 0.5*(mass*sum(sum(velocities.^2,2)))/num_particles;
        energy(energy_index,1) = energy(energy_index,1) - shifting_potential_term;
        energy(energy_index,3) = energy(energy_index,1) + energy(energy_index,2);
        steps = linspace(1, energy_index, energy_index);
        plot(ax1, steps, energy(1:energy_index,1), '-^', steps, energy(1:energy_index,2), '-v', steps, energy(1:energy_index,3), '-o');
        xlabel(ax1, 'Iterations');
        ylabel(ax1, 'Energy (\epsilon)');
        legend(ax1, 'Potential Energy', 'Kinetic Energy', 'Total Energy','Location','east');
    end
    
    if (mod(i, N_n) == 0)
         [neighbours_list, num_neighbours_list] = find_neighbours(num_particles, coordinates, length_cube, r_cutoff);
    end
    
    scatter3(ax2, coordinates(:, 1), coordinates(:, 2), coordinates(:, 3), 'filled');
    axis(ax2, [-0.5 10.5 -0.5 10.5 -0.5 10.5]   );
    pause(0.02);
end
toc;