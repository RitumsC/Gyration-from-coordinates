import get_gyr as gyr

molecule = gyr.GyrationMolecule("./molecules/water.xyz")
print("Data for: ", molecule.name)
print("Radiuses of gyration squared: ", molecule.Rx2, molecule.Ry2, molecule.Rz2, " [in angstroms**2]")
print("Gyration radius is: ", molecule.Rg, " [in angstroms]")
print("Relative shape anisotropy is: ", molecule.Rel_anis, " [dimensionless]")
print("Principal moments of inertia are: ", molecule.Ix2, molecule.Iy2, molecule.Iz2, " [in amu*angstroms^2 ]")
print("#####################################################################################\n\n")
molecule = gyr.GyrationMolecule("./molecules/methane.xyz")
print("Data for: ", molecule.name)
print("Radiuses of gyration squared: ", molecule.Rx2, molecule.Ry2, molecule.Rz2, " [in angstroms**2]")
print("Gyration radius is: ", molecule.Rg, " [in angstroms]")
print("Relative shape anisotropy is: ", molecule.Rel_anis, " [dimensionless]")
print("Principal moments of inertia are: ", molecule.Ix2, molecule.Iy2, molecule.Iz2, " [in amu*angstroms^2 ]")
print("#####################################################################################\n\n")
molecule = gyr.GyrationMolecule("./molecules/ether.xyz")
print("Data for: ", molecule.name)
print("Radiuses of gyration squared: ", molecule.Rx2, molecule.Ry2, molecule.Rz2, " [in angstroms**2]")
print("Gyration radius is: ", molecule.Rg, " [in angstroms]")
print("Relative shape anisotropy is: ", molecule.Rel_anis, " [dimensionless]")
print("Principal moments of inertia are: ", molecule.Ix2, molecule.Iy2, molecule.Iz2, " [in amu*angstroms^2 ]")
print("#####################################################################################\n\n")
molecule = gyr.GyrationMolecule("./molecules/isobutane.xyz")
print("Data for: ", molecule.name)
print("Radiuses of gyration squared: ", molecule.Rx2, molecule.Ry2, molecule.Rz2, " [in angstroms**2]")
print("Gyration radius is: ", molecule.Rg, " [in angstroms]")
print("Relative shape anisotropy is: ", molecule.Rel_anis, " [dimensionless]")
print("Principal moments of inertia are: ", molecule.Ix2, molecule.Iy2, molecule.Iz2, " [in amu*angstroms^2 ]")
print("#####################################################################################\n\n")
