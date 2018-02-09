function [x_new, v_new] = velocity_verlet(x, v, h, m, update_forces)
    a = (1/m) * update_forces(x);
    x_new = x + h*v + 0.5*h*h*a;
    v_new = v + 0.5*h*a;
    force = update_forces(x_new);
    v_new = v_new + 0.5*h*(1/m)*force;
end