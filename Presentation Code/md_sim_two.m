clear variables;
tic;
[num_particles, epsilon, sigma, r_cutoff, mass, density, kB, temperature, h, N_e, N_f, N_s, N_n, beta, gamma_one, gamma_two] = initialize_params_two();
length_cube = find_cube_length(num_particles, density);
diff = length_cube/((2*num_particles)^(1/3));
coordinates_one = initialize_cube_two(num_particles, length_cube-diff);
velocities_one = initialize_velocities(num_particles, mass, kB, temperature);
coordinates_two = coordinates_one;
velocities_two = velocities_one;


[neighbours_list_one, num_neighbours_list_one] = find_neighbours(num_particles, coordinates_one, length_cube, r_cutoff);
[forces_one, ~] = find_forces(num_particles, epsilon, sigma, coordinates_one, length_cube, neighbours_list_one, num_neighbours_list_one);

neighbours_list_two = neighbours_list_one;
num_neighbours_list_two = num_neighbours_list_one;
forces_two = forces_one;

energy_one = zeros(N_f/N_s, 3);
energy_two = energy_one;

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
    
    acceleration_one = forces_one/mass;
    coordinates_one = coordinates_one + h*velocities_one + 0.5*h^2*acceleration_one;
    velocities_one = velocities_one + h*(1-gamma_one)*acceleration_one;
    [forces_one, potential_energy_one] = find_forces(num_particles, epsilon, sigma, coordinates_one, length_cube, neighbours_list_one, num_neighbours_list_one);
    acceleration_one = forces_one/mass;
    velocities_one = velocities_one + h*gamma_one*acceleration_one;
    
    acceleration_two = forces_two/mass;
    coordinates_two = coordinates_two + h*velocities_two + 0.5*h^2*acceleration_two;
    velocities_two = velocities_two + h*(1-gamma_two)*acceleration_two;
    [forces_two, potential_energy_two] = find_forces(num_particles, epsilon, sigma, coordinates_two, length_cube, neighbours_list_two, num_neighbours_list_two);
    acceleration_two = forces_two/mass;
    velocities_two = velocities_two + h*gamma_two*acceleration_two;
    
    if (i <= N_e)
       scaling_factor_one = (3*kB*num_particles*temperature)/(mass*sum(sum(velocities_one.^2,2)));
       velocities_one = velocities_one * sqrt(scaling_factor_one);

       scaling_factor_two = (3*kB*num_particles*temperature)/(mass*sum(sum(velocities_two.^2,2)));
       velocities_two = velocities_two * sqrt(scaling_factor_two);
    elseif (mod(i, N_s) == 0)
        energy_index = (i-N_e)/N_s;
        energy_one(energy_index, 1) = potential_energy_one/num_particles;
        energy_one(energy_index, 2) = 0.5*(mass*sum(sum(velocities_one.^2,2)))/num_particles;
        energy_one(energy_index,1) = energy_one(energy_index,1) - shifting_potential_term;
        energy_one(energy_index,3) = energy_one(energy_index,1) + energy_one(energy_index,2);

        energy_two(energy_index, 1) = potential_energy_two/num_particles;
        energy_two(energy_index, 2) = 0.5*(mass*sum(sum(velocities_two.^2,2)))/num_particles;
        energy_two(energy_index,1) = energy_two(energy_index,1) - shifting_potential_term;
        energy_two(energy_index,3) = energy_two(energy_index,1) + energy_two(energy_index,2);
        
        steps = linspace(1, energy_index, energy_index);
        
        plot(ax1, steps, energy_one(1:energy_index,1), '-^', steps, energy_one(1:energy_index,2), '-v', steps, energy_one(1:energy_index,3), '-o');
        xlabel(ax1, 'Iterations');
        ylabel(ax1, 'Energy (\epsilon)');
        legend(ax1, 'Potential Energy', 'Kinetic Energy', 'Total Energy', 'Location', 'east');
        title(ax1, ['\beta = ', num2str(beta), ', \gamma = ', num2str(gamma_one)]); 
        
        plot(ax2, steps, energy_two(1:energy_index,1), '-^', steps, energy_two(1:energy_index,2), '-v', steps, energy_two(1:energy_index,3), '-o');
        xlabel(ax2, 'Iterations');
        ylabel(ax2, 'Energy (\epsilon)');
        legend(ax2, 'Potential Energy', 'Kinetic Energy', 'Total Energy', 'Location', 'east');
        title(ax2, ['\beta = ', num2str(beta), ', \gamma = ', num2str(gamma_two)]);
    end
    
    if (mod(i, N_n) == 0)
         [neighbours_list_one, num_neighbours_list_one] = find_neighbours(num_particles, coordinates_one, length_cube, r_cutoff);
         [neighbours_list_two, num_neighbours_list_two] = find_neighbours(num_particles, coordinates_two, length_cube, r_cutoff);
    end

    pause(0.02);
end
toc;