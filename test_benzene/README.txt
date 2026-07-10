Two sets are present, one for training and one for validation for a benzene system in a cubic box of 100 angstrom.
Here the format of the different input files is provided:

- energies : the list of energies in kcal/mol has to be provided;
- geometries : an extended .xyz where the first line is the number of atoms, the second line contains the three lattice vectors and the number of atomic species (kinds) in the system, the lines with the atomic coordinates where beyond the coordinates it is reported the atomic kind and the mass of the atom in atomic units;
- dipoles : the list of dipole moment values in atomic units (1x,1y,1z,2x,2y,...);
- forces : the list of gradients in atomic units (1x,1y,1z,2x,2y,..)

The files ending in *++ contain all the data including both training and validation.
