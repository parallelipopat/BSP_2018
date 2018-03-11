function [num_particles, epsilon, sigma, r_cutoff, mass, density, kB, temperature, h, N_e, N_f, N_s, beta, gamma] = initialize_params()
    num_particles = 864;
    epsilon = 1;
    sigma = 1;
    r_cutoff = 2.5;
    mass = 48;
    density = 38.74457642;
    kB = 1;
    temperature = 0.6694056515;
    h = 0.032;
    N_e = 10;
    N_f = 1000;
    N_s = 10;
    beta = 0;
    gamma = 0.2;
end