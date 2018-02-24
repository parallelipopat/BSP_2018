clear all; close all;

temperature = 94.4 / 119.5; % K
density = 0.8071593965; 
num_particles = 864;

N = 5;
h = 0.032;

[coordinates, length_cube] = initialize_cube(num_particles, density);
velocities = initialize_velocities(num_particles, temperature);
temp = zeros(N, 1);

beta = 0;
gamma = 0.5;

figure
hold on
for ii = 1 : N
    fprintf('%d\n', ii);
    %[coordinates, velocities] = newmark_beta(beta, gamma, coordinates, velocities, h, @update_acceleration, @d_acceleration);
    [coordinates, velocities] = velocity_verlet(coordinates, velocities, h, @update_acceleration);
    a = dot(velocities, velocities, 2);
    b = 16 * sum(a) / num_particles;
    temp(ii) = b;
    if (ii <= 2)
        velocities = (0.789958159 / b)^0.5 .* velocities;
    end
    plot(0:ii-1, temp(1:ii), '-o');
    pause(0.05);
end
%plot(1:N,temp,'-o');

function new_acceleration = update_acceleration(x)
    length_cube = 10.2294;
    LJ_force_reply = (1/48)*LJ_Force(x, length_cube);
    new_acceleration = LJ_force_reply(1);
end

function derivative_acceleration = d_acceleration(x)
    length_cube = 10.2294;
    LJ_force_reply = (1/48)*LJ_Force(x, length_cube);
    derivative_acceleration = LJ_force_reply(2);
end