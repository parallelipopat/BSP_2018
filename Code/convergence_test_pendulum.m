clear all;

% pendulum problem
l = 5; g = 10; T = 5;
theta = 0;
v_theta = 5;
y0 = [theta; v_theta];

f = @(y) [y(2); -g*(1/l)*sin(y(1))];

% reference solution with matlab function ode45
opts = odeset('AbsTol', 1e-12, 'RelTol', 1e-12);
[~, yref] = ode45(@(t,y) f(y), [0, T], y0, opts);
yref_final = yref(end,:).';

% main loop
N = 2.^(4:10); hsave = zeros(size(N));
err = zeros(size(N));
%beta = 0.5;
figure_num = 0;
%gamma = 0.5;
for gamma = 0:0.25:1
    for beta =0:0.25:1
        for ii = 1:length(N)
            h = T/N(ii); hsave(ii) = h;
            x = theta; v = v_theta;

            for jj = 1:N(ii)
                [x, v] = newmark_beta(beta, gamma, x, v, h, @update_acceleration, @d_acceleration);
                %[x, v] = newmark_beta_1d(beta, gamma, x, v, h, @update_acceleration, @d_acceleration);
            end

            err(ii) = norm([x; v] - yref_final);
        end

        figure_num = figure_num + 1;
        ax = subplot(5,5,figure_num);
        loglog(ax, hsave, err, '*-', 'linewidth', 2)
        fit = polyfit(log(hsave(end-3:end)), log(err(end-3:end)), 1);
        fprintf('The convergence order of the Newmark-Beta method with beta: %f and gamma: %f is %1.2f\n', beta, gamma, fit(1));
    end
end


function new_acceleration = update_acceleration(x)
    l = 5; g = 10;
    new_acceleration = -g*(1/l)*sin(x);
end

function derivative_acceleration = d_acceleration(x)
   l = 5; g = 10;
   derivative_acceleration = -g*(1/l)*cos(x);
end