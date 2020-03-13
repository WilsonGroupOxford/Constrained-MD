#ifndef FILE_IO_H
#define FILE_IO_H

#include "bond.h"
#include "constants.h"
#include <Eigen/Dense>
#include <memory>
#include <string>
#include <vector>

std::tuple<Eigen::MatrixXd, std::vector<std::string>> load_positions(const std::string& filename,
                                                                     const UnitType& unit_type,
                                                                     const int dimension);

std::vector<std::unique_ptr<Bond>> load_bonds(const std::string& filename,
                                              const Eigen::ArrayXXd& positions,
                                              const std::vector<std::string>& atom_names);

Eigen::VectorXd load_masses(const std::vector<std::string>& atom_types, const UnitType& unit_type);

void write_xyz(const std::string& filename, const Eigen::MatrixXd& positions,
               const Eigen::MatrixXd& velocities, const Eigen::MatrixXd& accelerations,
               const int step, const std::vector<std::string>& atom_types);

#endif // FILE_IO_H
