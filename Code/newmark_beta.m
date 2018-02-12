function [x_new, v_new] = newmark_beta(beta, x, v, h, update_acceleration)
    a = update_acceleration(x);
    x_new = x + h*v + 0.5*h*h*a;
    v_new = v + (1-beta)*h*a;
    a = update_acceleration(x_new);
    v_new = v_new + beta*h*a;
end