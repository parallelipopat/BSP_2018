clear all; close all;
[num_particles, epsilon, sigma, r_cutoff, mass, density, kB, temperature, h, N_e, N_f, N_s] = initialize_params();
length_cube = find_cube_length(num_particles, density);
diff = length_cube/((2*num_particles)^(1/3));
coordinates = initialize_cube(num_particles,length_cube-diff);
velocities = initialize_velocities(num_particles, mass, kB, temperature);
[neighbours_list,num_neighbours_list] = find_neighbours(num_particles, coordinates, length_cube, r_cutoff);

[forces, ~] = find_forces(num_particles, epsilon, sigma, coordinates, length_cube, neighbours_list, num_neighbours_list);
energy = zeros(N_f/N_s,3);

for i = 1:(N_e+N_f)
    fprintf('%d\n', i);
    acceleration = forces/mass;
    coordinates = coordinates + h*velocities + 0.5*h^2*acceleration;
    velocities = velocities + 0.5*h*acceleration;
    [forces, potential_energy] = find_forces(num_particles, epsilon, sigma, coordinates, length_cube, neighbours_list, num_neighbours_list);
    acceleration = forces/mass;
    velocities = velocities + 0.5*h*acceleration;
    %fprintf('%0.9f\n',velocities);
    if (i <= N_e)
       scaling_factor = (3*temperature*num_particles*kB)/(mass*sum(sum(velocities.^2,2)));
       velocities = velocities * sqrt(scaling_factor) ; 
    elseif (mod(i, N_s) == 0)
        energy_index = (i-N_e)/N_s;
        energy(energy_index, 1) = potential_energy/num_particles;
        energy(energy_index, 2) = 0.5*(mass*sum(sum(velocities.^2,2)))/num_particles;
        %fprintf('%0.9f\n',energy(energy_index,1));
        %fprintf('%0.9f\n\n',energy(energy_index,2));
    end       
end
energy(:,3) = energy(:,1) + energy(:,2);

steps = linspace(1,N_f/N_s,N_f/N_s);
figure;
plot(steps,energy(:,1),steps,energy(:,2),steps,energy(:,3));
xlabel('Steps'); %ylabel('Energy(\varepsilon)');
legend('Potential Energy', 'Kinetic Energy', 'Total Energy');