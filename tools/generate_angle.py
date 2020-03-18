import numpy as np
import sys
positions = dict()

if len(sys.argv) == 3:
    NO2_ANGLE = float(sys.argv[1])
    CO_ANGLE = float(sys.argv[2])
else:
    NO2_ANGLE = 0.0
    CO_ANGLE = np.pi / 2

with open("./CO_2_N2O_1.xyz") as fi:
    for line in fi:
        data = line.split()
        positions[data[0]] = np.array([float(item) for item in data[1:]])

bond_lengths = dict()

with open("./bonds.dat") as fi:
    fi.readline()
    for line in fi:
        atom_1, atom_2 = line.split()[1:3]
        separation = positions[atom_2] - positions[atom_1]
        distance = np.linalg.norm(separation)
        bond_lengths[f"{atom_1} {atom_2}"] = distance
        angle = np.arccos(separation[1] / distance)
        if separation[0] < 0:
            angle *= -1

desired_angles = {"Au%1 C%1": CO_ANGLE,
                  "C%1 O%1": CO_ANGLE,
                 "Au%1 C%2": -np.pi/2,
                  "C%2 O%2": -np.pi/2,
                  "Au%1 N%2": NO2_ANGLE,
                  "N%2 N%1": NO2_ANGLE,
                  "N%1 O%3": NO2_ANGLE}

new_positions = {"Au%1": np.array([0.0, 0.0, 0.0])}
for key, val in desired_angles.items():
    atom_1, atom_2 = key.split()
    new_positions[atom_2] = new_positions[atom_1] + bond_lengths[key] * np.array([np.sin(val), np.cos(val), 0.0])

with open("./angle.xyz", "w") as fi:
    fi.write(f"{len(new_positions)}\n")
    fi.write(f"Shifted by {CO_ANGLE}, and {NO2_ANGLE}\n")
    for key, val in new_positions.items():
        fi.write(f"{key} {val[0]} {val[1]} {val[2]}\n")

