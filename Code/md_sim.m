clear all; close all;

[num_particles, epsilon, sigma, mass, density, kB, temperature, h, N] = initialise_params();

length_cube = find_cube_length(num_particles, density);
%coordinates = initialize_cube(num_particles, length_cube);
coordinates = initialize_position(864,3,4,[6,6,6],5.45*[1,1,1]);
coordinates_array = zeros(num_particles, 3, N+1);
velocities_array = zeros(num_particles, 3, N+1);
acceleration_array = zeros(num_particles, 3, N+1);
velocities = initialize_velocities(num_particles, temperature);
temp_array = zeros(N, 1);

coordinates_array(:,:,1) = coordinates;
velocities_array(:,:,1) = coordinates;
acceleration_array(:,:,1) = update_acceleration(coordinates);

figure
hold on
for ii = 1 : N
    fprintf('%d ', ii);
    %[coordinates, velocities] = newmark_beta(beta, gamma, coordinates, velocities, h, @update_acceleration, @d_acceleration);
    [coordinates, velocities] = velocity_verlet(coordinates, velocities, h, @update_acceleration);
    b = mass * sum(dot(velocities, velocities, 2)) / (3*num_particles*kB);
       
    if (ii <= 2)
        velocities = (temperature / b)^0.5 .* velocities;
    end
    b = mass * sum(dot(velocities, velocities, 2)) / (3*num_particles*kB);
    fprintf('%.9f, ', b);
    temp_array(ii) = b;
    coordinates_array(:,:,ii+1) = coordinates;
    velocities_array(:,:,ii+1) = velocities;
    acceleration_array(:,:,ii+1) = update_acceleration(coordinates);
    ax1 = subplot(1,2,1);
    ax2 = subplot(1,2,2);
    plot(ax1, 0:ii-1, temp_array(1:ii), '-o');
    scatter3(ax2, coordinates(:, 1), coordinates(:, 2), coordinates(:, 3), 'filled');
    pause(0.05);
end
%plot(1:N,temp,'-o');

function new_acceleration = update_acceleration(x)
    [num_particles, ~, ~, mass, density, ~, ~, ~, ~] = initialise_params();
    length_cube = find_cube_length(num_particles, density);
    new_acceleration = (1/mass)*LJ_Force(x, length_cube);
end

% function derivative_acceleration = d_acceleration(x)
%     length_cube = 0.0071;
%     LJ_force_reply = (1/48)*LJ_Force(x, length_cube);
%     derivative_acceleration = LJ_force_reply(2);
% end