function vector = periodic_boundary_correction_1d(vector, length_cube)
    half_box_length = length_cube / 2;
    
    if (vector > half_box_length)
            vector = vector - length_cube;
    elseif (vector < half_box_length)
            vector = vector + length_cube;
    end
end