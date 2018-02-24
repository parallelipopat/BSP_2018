function [x_new, v_new] = newmark_beta_1d(beta, gamma, x_i, v_i, h, update_acceleration, d_acceleration)
    a_i = update_acceleration(x_i);
    v_new = v_i + h*(1-gamma)*a_i; 
    
    f = @(z, v_i, a_i, h) v_i + 0.5*h*((1-2*beta)*a_i + 2*beta*update_acceleration(z));
    Df = @(z, h) 0.5*beta*h*h*z;

    F = @(x_j, x_i, h) x_j - x_i - h*f(x_j, v_i, a_i, h);
    DF = @(x_j, h) eye(size(x_j, 1)) - Df(d_acceleration(x_j), h);
    N = @(x_i, h) x_i+h - (DF(x_i+h, h))^-1*F(x_i+h, x_i, h);
    
    x_new = N(x_i, h); 
    a_new = update_acceleration(x_new); 
    v_new = v_new + h*gamma*a_new;
end