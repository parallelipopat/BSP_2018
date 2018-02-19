clear all; close all;

temperature = 94.4 / 119.5; % K
density = 0.8071593965; 
num_particles = 864;

N = 100;
h = 0.032;

[coordinates, length_cube] = initialize_cube(num_particles, density);
velocities = initialize_velocities(num_particles, temperature);
temp = zeros(N, 1);

for ii = 1 : N
    fprintf('%d\n', ii);
    [coordinates, velocities] = velocity_verlet(coordinates, velocities, h, @update_acceleration);
    a = dot(velocities, velocities, 2);
    b = 16 * sum(a) / num_particles;
    temp(ii) = b;
    if (ii <= 50)
        velocities = (0.789958159 / b)^0.5 .* velocities;
    end
end
plot(temp);

function new_acceleration = update_acceleration(x)
    length_cube = 10.2294;
    new_acceleration = LJ_Force(x, length_cube);
end