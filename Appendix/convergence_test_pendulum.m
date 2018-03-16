% Code adapted from A. Paganini, Numerical Solutions of DEs I.

% Pendulum problem initial configuration
l = 5; g = 10; T = 5;
theta = 0;
v_theta = 5;
y0 = [theta; v_theta];

f = @(y) [y(2); -g*(1/l)*sin(y(1))];

% Reference solution with MATLAB function ode45
opts = odeset('AbsTol', 1e-12, 'RelTol', 1e-12);
[~, yref] = ode45(@(t,y) f(y), [0, T], y0, opts);
yref_final = yref(end,:).';

N = 2.^(4:10); hsave = zeros(size(N));
err = zeros(size(N));
figure_num = 0;
beta = 1.00; % we vary this value as preferred

% For each value of gamma, repeat the pendulum problem for different step sizes h
for gamma =0:0.25:1
    for i = 1:length(N)
        h = T/N(i); hsave(i) = h;
        x = theta; v = v_theta;

        for jj = 1:N(i)
            [x, v] = newmark_beta(beta, gamma, x, v, h, ...
                                  @update_acceleration, @d_acceleration);
         end
         
        err(i) = norm([x; v] - yref_final);
    end

    % Plot the graph of convergence
    loglog(hsave, err, '*-', 'linewidth', 2)
    ylim([10^-5 10^1]);
    xlabel('log(h)');
    ylabel('log(err)');
    title('\beta = 1.00');
    hold on;
    % Find (asymptotic) slope of the graph plotted
    fit = polyfit(log(hsave(end-3:end)), log(err(end-3:end)), 1);
    fprintf('The convergence order of the Newmark-Beta method with \beta: %f ...
             and \gamma: %f is %1.2f\n', beta, gamma, fit(1));
end
hold off;
legend('0.00','0.25','0.50','0.75','1.00');
title(legend,'\gamma values');
legend('Location','northwest');




% Define the acceleration function and its derivate (for Newton-Raphson iteration)
function new_acceleration = update_acceleration(x)
    l = 5; g = 10;
    new_acceleration = -g*(1/l)*sin(x);
end

function derivative_acceleration = d_acceleration(x)
   l = 5; g = 10;
   derivative_acceleration = -g*(1/l)*cos(x);
end