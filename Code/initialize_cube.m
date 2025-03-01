function r = initialize_position(num_particles, length_cube) % FCC crystal
    one_side = (2*num_particles)^(1/3); 
    r = zeros(num_particles,3);
    i = 1;
    particle_spacing = length_cube / (one_side - 1);
    for x = 0:particle_spacing:length_cube
        for y = 0:particle_spacing:length_cube
            if (mod(x+y,2*particle_spacing) == 0)
               for z = 0:2*particle_spacing:length_cube
                  r(i,:) = [x,y,z];
                  i = i + 1;
               end
            else
                for z = particle_spacing:2*particle_spacing:length_cube
                  r(i,:) = [x,y,z];
                  i = i + 1;
               end
            end
        end
    end
end

% function coordinates = initialize_cube(num_particles, length_cube) 
%     seed = 7;
%     rng(seed);
%     
%     a = -length_cube/2;
%     b = -a;
%     %a = 0; b = length_cube;
%     coordinates = (b-a).*rand(num_particles, 3)+ a;
%     %scatter3(coordinates(:, 1), coordinates(:, 2), coordinates(:, 3), 'filled');
% end