function vector = periodic_boundary_correction(vector, length_cube)
    half_box_length = length_cube / 2;
    
    for dimension = 1:3
        if (vector(dimension) > half_box_length)
            vector(dimension) = vector(dimension) - length_cube;
        elseif (vector(dimension) < -half_box_length)
            vector(dimension) = vector(dimension) + length_cube;
        end
    end
end