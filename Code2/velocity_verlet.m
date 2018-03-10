function [x_new, v_new] = velocity_verlet(x, v, h, update_acceleration)
    a = update_acceleration(x);
    x_new = x + h*v + 0.5*h*h*a;
    v_new = v + 0.5*h*a;
    a = update_acceleration(x_new);
    v_new = v_new + 0.5*h*a;
end