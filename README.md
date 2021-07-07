# Gyration-from-coordinates
Calculates gyration, relative shape anisotropy and moments of inertia from coordinate files(e.g. xyz).

Module to calculate gyaration and inertial information from coordinates.

Currently only .xyz file input is supported.

To use:
1) import the module in your python script with (script must be in same directory as this file)
import get_gyr as gyr
2) read .xyz file
molecule = gyr.GyrationMolecule("./molecules/water.xyz")
3) call appropriate variable, e.g.
print("Relative shape anisotropy is: ", molecule.Rel_anis)

for detailed information see test_gyr.py script
