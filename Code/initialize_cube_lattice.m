function coordinates = initialize_cube_lattice(num_mol, density)
    coordinates = zeros(num_mol, 3);
    num_side_molecules = ceil((num_mol)^(1/3));
    length = (num_mol/density)^(1/3);
    index = [0; 0; 0];
    
    for mol = 1 : num_mol
       coordinates(mol,:) = (index + [0.5; 0.5; 0.5])*(length/num_side_molecules);
       
       index(3) = index(3) + 1;
       if (index(3) == num_side_molecules)
           index(3) = 0;
           index(2) = index(2) + 1;
           if (index(2) == num_side_molecules)
               index(2) = 0;
               index(1) = index(1) + 1;
           end
       end 
    end
    
    figure;
    scatter3(coordinates(:, 1), coordinates(:, 2), coordinates(:, 3), 'filled');
    view(-77,17);
end