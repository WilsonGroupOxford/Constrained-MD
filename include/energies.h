#pragma once

#include "energies.h"

#include "bond.h"

#include <Eigen/Dense>
#include <memory>
#include <vector>

double calculate_kinetic_energy(const Eigen::MatrixXd& velocities, const Eigen::MatrixXd& masses);

double calculate_potential_energy(const std::vector<std::unique_ptr<Bond>>& bonds, const Eigen::MatrixXd& positions);
