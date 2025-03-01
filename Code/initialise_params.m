function [num_particles, epsilon, sigma, mass, density, kB, temperature, h, N] = initialise_params()
    num_particles = 864;
    epsilon = 1;
    sigma = 1;
    mass = 48;
    density = 38.74457642; %1.374; % g cm^-3 0.8071593965;
    kB = 1; %J K^-1
    temperature = 94.4/119.509;
    h = 0.032; %s
    N = 20;

%     num_particles = 864;
%     epsilon = 1.65 * 10^(-21);
%     sigma = 3.4 * 10^(-10);
%     mass = 39.95 * 1.6747 * 10^(-27); %g
%     density = 1.374; % g cm^-3 0.8071593965;
%     kB = 1.38064852 * 10^(-23); %J K^-1
%     temperature = 94.4;
%     h = 10^(14); %s
%     N = 20;

end