function coordinates = initialize_cube(num_particles, length_cube) % FCC crystal
    one_side = (2*num_particles)^(1/3); 
    coordinates = zeros(3*num_particles,1);
    i = 1;
    particle_spacing = length_cube / (one_side - 1);
    for x = 0:particle_spacing:length_cube
        for y = 0:particle_spacing:length_cube
            if (mod(x+y,2*particle_spacing) == 0)
               for z = 0:2*particle_spacing:length_cube
                  coordinates(i:i+2) = [x,y,z]';
                  i = i + 3;
               end
            else
                for z = particle_spacing:2*particle_spacing:length_cube
                  coordinates(i:i+2) = [x,y,z]';
                  i = i + 3;
               end
            end
        end
    end
end
