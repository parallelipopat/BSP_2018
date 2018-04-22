function [num_particles, epsilon, sigma, r_cutoff, mass, density, kB, temperature, h, N_e, N_f, N_s, N_n, beta, gamma_one, gamma_two] = initialize_params_two()
    num_particles = 864;
    epsilon = 1;
    sigma = 1;
    r_cutoff = 2.5;
    mass = 48;
    density = 38.74457642;
    kB = 1;
    temperature = 1;
    h = 0.032;
    N_e = 0;
    N_f = 100;
    N_s = 5;
    N_n = 15;
    beta = 0.0;
    gamma_one = 0.20; % between 0.0 and 0.5
    gamma_two = 0.70; % between 0.5 and 1.0
end