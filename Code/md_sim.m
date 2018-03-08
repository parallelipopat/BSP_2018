clear all; close all;

[num_particles, epsilon, sigma, mass, density, kB, temperature, h] = initialise_params()
N = 20;

main_array = zeros(num_particles,3,N);
length_cube = find_cube_length(num_particles, density);
coordinates = initialize_cube(num_particles, length_cube);
coordinates_array = zeros(num_particles, 3, N);
velocities = initialize_velocities(num_particles, temperature);
temp_array = zeros(N, 1);

figure
hold on
for ii = 1 : N
    fprintf('%d\n', ii);
    %[coordinates, velocities] = newmark_beta(beta, gamma, coordinates, velocities, h, @update_acceleration, @d_acceleration);
    a = dot(velocities, velocities, 2);
    b = mass * sum(a) / (3*num_particles*kB);
    temp_array(ii) = b;
    main_array(:,:,ii) = coordinates;
    [coordinates, velocities] = velocity_verlet(coordinates, velocities, h, @update_acceleration);
    coordinates_array(:,:,ii) = coordinates;
    if (ii <= 5)
        velocities = (temperature / b)^0.5 .* velocities;
    end
    plot(0:ii-1, temp_array(1:ii), '-o');
    pause(0.05);
end
%plot(1:N,temp,'-o');

function new_acceleration = update_acceleration(x)
    [num_particles, epsilon, sigma, mass, density, kb, temperature, h] = initialise_params();
    length_cube = find_cube_length(num_particles, density);
    LJ_force_reply = (1/mass)*LJ_Force(x, length_cube);
    new_acceleration = LJ_force_reply(1);
end

function derivative_acceleration = d_acceleration(x)
    length_cube = 0.0071;
    LJ_force_reply = (1/48)*LJ_Force(x, length_cube);
    derivative_acceleration = LJ_force_reply(2);
end