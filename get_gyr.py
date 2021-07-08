"""Module to calculate gyaration and inertial information from coordinates.

Currently only .xyz file input is supported.

To use:
1) import the module in your python script with (script must be in same directory as this file)
import get_gyr as gyr
2) read .xyz file
molecule = gyr.GyrationMolecule("./molecules/water.xyz")
3) call appropriate variable, e.g.
print("Relative shape anisotropy is: ", molecule.Rel_anis)

for detailed information see test_gyr.py script
"""

import numpy as np
import re


class GyrationMolecule:

    def __init__(self, fname, flag_xyz=1, flag_periodic=0):
        if flag_xyz:
            self.name = fname.split('.xyz')[0]

            # Gyration parameters
            self.__coords = self.__xyz_to_coords(fname)
            self.Rx2, self.Ry2, self.Rz2, = self.__calc_evals_gyr(self.__coords)
            self.Rg = np.sqrt(self.Rx2 + self.Ry2 + self.Rz2)
            self.Rel_anis = (1.5 * (self.Rx2**2 + self.Ry2**2 + self.Rz2**2)
                             / (self.Rx2 + self.Ry2 + self.Rz2)**2 - 0.5)

            # Inertia parametrs
            self.Ix2, self.Iy2, self.Iz2, = self.__calc_evals_inertia(self.__coords)
        else:
            print("Not yet implemented")
            quit()

    __massDictionary = {
        "H": 1.008,
        "He": 4.003,
        "Li": 6.941,
        "Be": 9.012,
        "B": 10.811,
        "C": 12.011,
        "N": 14.007,
        "O": 15.999,
        "F": 18.998,
        "Ne": 20.18,
        "Na": 22.99,
        "Mg": 24.305,
        "Al": 26.982,
        "Si": 28.086,
        "P": 30.974,
        "S": 32.065,
        "Cl": 35.453,
        "Ar": 39.948,
        "K": 39.098,
        "Ca": 40.078,
        "Sc": 44.956,
        "Ti": 47.867,
        "V": 50.942,
        "Cr": 51.996,
        "Mn": 54.938,
        "Fe": 55.845,
        "Co": 58.933,
        "Ni": 58.693,
        "Cu": 63.546,
        "Zn": 65.39,
        "Ga": 69.723,
        "Ge": 72.64,
        "As": 74.922,
        "Se": 78.96,
        "Br": 79.904,
        "Kr": 83.8,
        "Rb": 85.468,
        "Sr": 87.62,
        "Y": 88.906,
        "Zr": 91.224,
        "Nb": 92.906,
        "Mo": 95.94,
        "Tc": 98,
        "Ru": 101.07,
        "Rh": 102.906,
        "Pd": 106.42,
        "Ag": 107.868,
        "Cd": 112.411,
        "In": 114.818,
        "Sn": 118.71,
        "Sb": 121.76,
        "Te": 127.6,
        "I": 126.905,
        "Xe": 131.293,
        "Cs": 132.906,
        "Ba": 137.327,
        "La": 138.906,
        "Ce": 140.116,
        "Pr": 140.908,
        "Nd": 144.24,
        "Pm": 145,
        "Sm": 150.36,
        "Eu": 151.964,
        "Gd": 157.25,
        "Tb": 158.925,
        "Dy": 162.5,
        "Ho": 164.93,
        "Er": 167.259,
        "Tm": 168.934,
        "Yb": 173.04,
        "Lu": 174.967,
        "Hf": 178.49,
        "Ta": 180.948,
        "W": 183.84,
        "Re": 186.207,
        "Os": 190.23,
        "Ir": 192.217,
        "Pt": 195.078,
        "Au": 196.967,
        "Hg": 200.59,
        "Tl": 204.383,
        "Pb": 207.2,
        "Bi": 208.98,
        "Po": 209,
        "At": 210,
        "Rn": 222,
        "Fr": 223,
        "Ra": 226,
        "Ac": 227,
        "Th": 232.038,
        "Pa": 231.036,
        "U": 238.029,
        "Np": 237,
        "Pu": 244,
        "Am": 243,
        "Cm": 247,
        "Bk": 247,
        "Cf": 251,
        "Es": 252,
        "Fm": 257,
        "Md": 258,
        "No": 259,
        "Lr": 262,
        "Rf": 261,
        "Db": 262,
        "Sg": 266,
        "Bh": 264,
        "Hs": 277,
        "Mt": 268
    }

    @staticmethod
    def __xyz_to_coords(xyz_file: str):
        """ Returns atom mass and coordinates from xyz file

        Arguments:
            xyz_file: location of xyz file to convert to coordinates
        Returns:
            atom_coordinates: [[mass1, ,mass2,...],[x1, x2...],[y1, y2...],[z1, z2...]]
        """
        with open(xyz_file) as file:
            atom_coordinates = [[], [], [], []]
            i = 0
            for line in file:
                if i >= 2:
                    atom, x, y, z = re.split(' +|\t|\n', line)[0:-1]
                    atom_coordinates[0].append(float(GyrationMolecule.__massDictionary[atom]))
                    atom_coordinates[1].append(float(x))
                    atom_coordinates[2].append(float(y))
                    atom_coordinates[3].append(float(z))
                i += 1

        return atom_coordinates

    @staticmethod
    def __calc_com_gyr(atom_coordinates, periodic: int=0) -> (float, float, float):
        """Returns center of points given atom coordinates under specified conditions.

        Arguments:
            atom_coordinates: [[mass1, ,mass2,...],[x1, x2...],[y1, y2...],[z1, z2...]]
            periodic: 1 if coordinates given in periodic system, else 0
        Returns:
            x_com, y_com, z_com: center of mass for each coordinate
        """
        if periodic:
            print("Gyration calculations for periodic conditions not yet implemented!")
            quit()
        else:
            a_avg = [0.0, 0.0, 0.0]
            for i in range(1, 4):
                for coord in atom_coordinates[i]:
                    a_avg[i-1] += coord
                a_avg[i-1] = a_avg[i-1] / len(atom_coordinates[1])
            x_com, y_com, z_com = a_avg
        return x_com, y_com, z_com

    @staticmethod
    def __calc_com_inertia(atom_coordinates, periodic: int=0) -> (float, float, float):
        """Returns center of mass given atom coordinates under specified conditions.

        Arguments:
            atom_coordinates: [[mass1, ,mass2,...],[x1, x2...],[y1, y2...],[z1, z2...]]
            periodic: 1 if coordinates given in periodic system, else 0
        Returns:
            x_com, y_com, z_com: center of mass for each coordinate
        """
        if periodic:
            print("Inertia calculations for periodic conditions not yet implemented!")
            quit()
        else:
            a_avg = [0.0, 0.0, 0.0]
            total_mass = 0.0
            for i in range(1, 4):
                atom_index = 0
                for coord in atom_coordinates[i]:
                    a_avg[i-1] += coord*atom_coordinates[0][atom_index]
                    total_mass += atom_coordinates[0][atom_index]
                    atom_index += 1
                a_avg[i-1] = a_avg[i-1] / total_mass
            x_com, y_com, z_com = a_avg
        return x_com, y_com, z_com

    @staticmethod
    def __calc_evals_gyr(atom_coordinates, periodic:int=0) -> (float, float, float):
        """Returns eigenvalues squared for gyration tensor.
            For details refer to: https://en.wikipedia.org/wiki/Gyration_tensor

        Arguments:
            atom_coordinates: [[mass1, ,mass2,...],[x1, x2...],[y1, y2...],[z1, z2...]]
            periodic: 1 if coordinates given in periodic system, else 0
        Returns:
            Rx2, Ry2, Rz2: squared eigenvalues for all axis [AA**2]
        """
        rx, ry, rz = 0.0, 0.0, 0.0

        if periodic:
            print("Gyration for periodic systems not yet implemented!")
            quit()
        else:
            # Calculate centres of mass and shift coordinates by it
            coms = GyrationMolecule.__calc_com_gyr(atom_coordinates, periodic)
            for i in range(1, 4):
                for coord in atom_coordinates[i]:
                    coord -= coms[i-1]

            # Calculate individual matrix components
            def calc_elem(axis1, axis2):
                axis_dot = np.dot(axis1, axis2)
                return axis_dot / len(axis1)

            s_xx = calc_elem(atom_coordinates[1], atom_coordinates[1])
            s_xy = calc_elem(atom_coordinates[1], atom_coordinates[2])
            s_xz = calc_elem(atom_coordinates[1], atom_coordinates[3])
            s_yy = calc_elem(atom_coordinates[2], atom_coordinates[2])
            s_yz = calc_elem(atom_coordinates[2], atom_coordinates[3])
            s_zz = calc_elem(atom_coordinates[3], atom_coordinates[3])

            # Form matrix
            r_matrix = [[s_xx, s_xy, s_xz],
                        [s_xy, s_yy, s_yz],
                        [s_xz, s_yz, s_zz]]

            # Diagonalise the matrix and return eigenvalues
            rx, ry, rz = np.linalg.eigvals(r_matrix)

        return rx * rx, ry * ry, rz * rz

    @staticmethod
    def __calc_evals_inertia(atom_coordinates, periodic:int=0) -> (float, float, float):
        """Returns eigenvectors for inertia tensor.

        For details refer to: https://en.wikipedia.org/wiki/Moment_of_inertia#Inertia_tensor

        Arguments:
            atom_coordinates: [[mass1, ,mass2,...],[x1, x2...],[y1, y2...],[z1, z2...]]
            periodic: 1 if coordinates given in periodic system, else 0
        Returns:
            Ix, Iy, Iz: eigenvalues for inertia tensor's all axis [amu*AA^2]
        """
        rx, ry, rz = 0.0, 0.0, 0.0

        if periodic:
            print("Inertia in periodic systems - Not yet implemented!")
            quit()
        else:
            # Calculate centres of mass and shift coordinates by it
            coms = GyrationMolecule.__calc_com_inertia(atom_coordinates, periodic)
            for i in range(1, 4):
                for coord in atom_coordinates[i]:
                    coord -= coms[i-1]

            # Calculate individual matrix components
            def calc_diag_elem(axis2, axis3):
                i_ii = 0
                for i in range(0, len(axis2)):
                    i_ii += atom_coordinates[0][i]*(axis2[i]**2 + axis3[i]**2)
                return i_ii

            def calc_prod_inert(axis1, axis2):
                i_ij = 0
                for i in range(0, len(axis2)):
                    i_ij -= atom_coordinates[0][i] * axis1[i]*axis2[i]
                return i_ij

            i_xx = calc_diag_elem(atom_coordinates[2], atom_coordinates[3])
            i_xy = calc_prod_inert(atom_coordinates[1], atom_coordinates[2])
            i_xz = calc_prod_inert(atom_coordinates[1], atom_coordinates[3])
            i_yy = calc_diag_elem(atom_coordinates[1], atom_coordinates[3])
            i_yz = calc_prod_inert(atom_coordinates[2], atom_coordinates[3])
            i_zz = calc_diag_elem(atom_coordinates[1], atom_coordinates[2])

            # Form matrix
            i_matrix = [[i_xx, i_xy, i_xz],
                        [i_xy, i_yy, i_yz],
                        [i_xz, i_yz, i_zz]]

            # Diagonalise the matrix and return eigenvalues
            ix, iy, iz = np.linalg.eigvals(i_matrix)

        return ix, iy, iz

