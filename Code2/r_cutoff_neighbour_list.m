[num_particles, epsilon, sigma, r_cutoff, mass, density, kB, temperature, h, N_e, N_f, N_s, beta, gamma] = initialize_params();
length_cube = find_cube_length(num_particles, density);
coordinates = initialize_cube(num_particles, length_cube); 
distance_array = zeros(0.5*num_particles*(num_particles-1),1);
k = 1;
for i = 1:num_particles-1
        for j = i+1:num_particles
            diff_r = coordinates(i,:) - coordinates(j, :);
            diff_r = diff_r - length_cube*round(diff_r/length_cube);
            norm_r = norm(diff_r);
            distance_array(k) = round(norm_r,3);
            k = k+1;
        end
end
% dist_r = 2.00:0.01:3.00;
% dist_r_6 = dist_r.^(-6);
% potential = 4*(dist_r_6).*((dist_r_6) - 1);

% plot(dist_r, potential);
histogram(distance_array);