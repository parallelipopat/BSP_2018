clear all;

% pendulum problem
l = 5; g = 10; T = 5; m =1;
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
err = zeros(size(N)); err_EE = zeros(size(N));
for ii = 1:length(N)
    h = T/N(ii); hsave(ii) = h;
    x = theta; v = v_theta; ysol_EE = y0; 
    
    for jj = 1:N(ii)
        [x, v] = velocity_verlet(x, v, h, m, @update_forces);
        ysol_EE = ysol_EE + h*f(ysol_EE);
        %keyboard
    end
    
    err(ii) = norm([x; v] - yref_final);
    %keyboard
    err_EE(ii) = norm(ysol_EE - yref_final);
end

figure;
loglog(hsave, err_EE, '*-', 'linewidth', 2)
hold on
loglog(hsave, err, '*-', 'linewidth', 2)
legend('EE', 'VV');
fit = polyfit(log(hsave), log(err), 1);
fit_EE = polyfit(log(hsave), log(err_EE), 1);
fprintf('the convergence order of Velocity Verlet is %1.2f\n', fit(1));
fprintf('the convergence order of Explicit Euler is %1.2f\n', fit_EE(1));

function new_force = update_forces(x)
    l = 5; g = 10;
    new_force = -g*(1/l)*sin(x);
end