function [num_particles, epsilon, sigma, r_cutoff, mass, density, kB, temperature, h, N_e, N_f, N_s] = initialize_params()
    num_particles = 864;
    epsilon = 1;
    sigma = 1;
    r_cutoff = 2.5;
    mass = 48;
    density = 38.74457642; %1.374; % g cm^-3 0.8071593965;
    kB = 1; %J K^-1
    temperature = 80.0/119.509;
    h = 0.032; %s
    N_e = 5;
    N_f = 20;
    N_s = 1;
    
%     num_particles = 864;
%     epsilon = 1.65 * 10^(-21);
%     sigma = 3.4 * 10^(-10);
%     r_cutoff = 2.5 * sigma;
%     mass = 39.95 * 1.6747 * 10^(-27); %g
%     density = 1.374; % g cm^-3 0.8071593965;
%     kB = 1.38064852 * 10^(-23); %J K^-1
%     temperature = 94.4;
%     h = 10^(14); %s
%     N_e = 5;
%     N_p = 20;
%     N_s = 2;
end