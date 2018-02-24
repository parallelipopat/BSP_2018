function [x_new, v_new] = newmark_beta(beta, gamma, x_i, v_i, h, update_acceleration, d_acceleration)
    a_i = update_acceleration(x_i);
    y_i = [v_i; x_i];

    
    f = @(z, v_i, a_i, h) [(1-gamma)*a_i + gamma*update_acceleration(z); v_i + 0.5*h*((1-2*beta)*a_i + 2*beta*update_acceleration(z))];
    Df = @(z, h) [zeros(size(z)) h*gamma*z; zeros(size(z)) 0.5*beta*h*h*z];

    F = @(y_j, y_i, h) y_j - y_i - h*f(y_j(2), v_i, a_i, h);
    DF = @(y_j, h) eye(size(y_j, 1)) - Df(d_acceleration(y_j(2)), h);
    N = @(y_i, h) y_i+h - (DF(y_i+h, h))^-1*F(y_i+h,y_i,h);
    
    y_j = N(y_i, h); 
    x_new = y_j(2);
    v_new = y_j(1);

end