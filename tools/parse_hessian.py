#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 11 13:53:30 2020

@author: newc4592
"""

import numpy as np

def load_coordinates(filename: str):
    coordinates = []
    atom_names = []
    with open(filename) as fi:
        fi.readline()
        fi.readline()
        for line in fi:
            if not line:
                break
            line = line.split()
            atom_names.append(line[0])
            coordinates.append(np.array([float(item) for item in line[1:]]))
    return atom_names, coordinates

def parse_hessian(filename: str):
    data = np.genfromtxt(filename, dtype=float).ravel()

    # Gaussian stores the Hessian as a lower triangular matrix.
    # We can solve 3X(3X+1)/2 = N using the quadratic formula:
    degrees_of_freedom = (-1 + np.sqrt(1 + 8 * np.shape(data)[0])) / 2
    degrees_of_freedom = int(degrees_of_freedom)
    hessian = np.empty([degrees_of_freedom, degrees_of_freedom], dtype=float)
    tril_indices = np.tril_indices(degrees_of_freedom)

    for count in range(tril_indices[0].shape[0]):
        i, j = tril_indices[0][count], tril_indices[1][count]
        hessian[i, j] = data[count]
        hessian[j, i] = data[count]
    return hessian


def hessian_to_force_constants(hessian, coordinates, bonds):
    degrees_of_freedom = hessian.shape[0]

    # Calculate the stretches for all pairs of molecules

    force_constants = np.zeros([degrees_of_freedom // 3, degrees_of_freedom // 3])
    for bond in bonds:
        between_atom_vec = coordinates[bond[1]] - coordinates[bond[0]]
        between_atom_vec /= np.linalg.norm(between_atom_vec)
        d_vec = np.zeros(degrees_of_freedom)
        for other_atom in bond[1:]:
            d_vec[3*(other_atom): 3 * (other_atom + 1)] = between_atom_vec
        force_constant = np.matmul(np.transpose(d_vec), np.matmul(hessian, d_vec))
        force_constants[bond[0], bond[1]] = force_constant
        force_constants[bond[1], bond[0]] = force_constant

    return force_constants


def coordinates_to_distances(coordinates):
    atoms = len(coordinates)

    equilibrium_distances = np.empty([atoms, atoms])
    for atom in range(atoms):
        equilibrium_distances[atom, atom] = 0.0
        for other_atom in range(atom):
            distance = np.linalg.norm(coordinates[atom] - coordinates[other_atom])
            equilibrium_distances[atom, other_atom] = distance
            equilibrium_distances[other_atom, atom] = distance
    return equilibrium_distances


def write_bond_file(bonds, force_constants, distances, atom_names):
    with open("./bonds.dat", "w") as fi:
        fi.write("arbitrary");
        for bond in bonds:
            i = atom_names.index(bond[0])
            j = atom_names.index(bond[1])
            fi.write(f"{bond[0]} {bond[1]} {force_constants[i, j]:.5f} {distances[i, j]:.5f}\n")


if __name__ == "__main__":
    BONDED_PAIRS = [("Au", "C1", "O1"),
                    ("C1", "O1"),
                    ("Au", "C2", "O2"),
                    ("C2", "O2"),
                    ("Au", "N2", "N1", "O3"),
                    ("N2", "N1", "O3"),
                    ("N1", "O3")]
    ATOM_NAMES, COORDINATES = load_coordinates("CO_2_N2O_1.xyz")
    HESSIAN = parse_hessian("CO_2_N2O_1.dat")
    BOND_IDS = [tuple(ATOM_NAMES.index(NAME) for NAME in BOND) for BOND in BONDED_PAIRS]
    FORCE_CONSTANTS = hessian_to_force_constants(HESSIAN, COORDINATES, BOND_IDS)
    DISTANCES = coordinates_to_distances(COORDINATES)
    write_bond_file(BONDED_PAIRS, FORCE_CONSTANTS, DISTANCES, ATOM_NAMES)