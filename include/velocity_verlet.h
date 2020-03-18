#ifndef VELOCITY_VERLET_H
#define VELOCITY_VERLET_H

#include <Eigen/Dense>
#include <memory>
#include <vector>

#include "bond.h"
#include "constants.h"


Eigen::MatrixXd calculate_forces(const std::vector<std::unique_ptr<Bond>>& bonds,
                                 const Eigen::MatrixXd& positions);

Eigen::MatrixXd velocity_verlet(Eigen::MatrixXd& positions, Eigen::MatrixXd& velocities,
                                Eigen::MatrixXd& accelerations, const Eigen::VectorXd& masses,
                                const double dt, const std::vector<std::unique_ptr<Bond>>& bonds);

void excite_bond(const std::unique_ptr<Bond>& bond, Eigen::MatrixXd& positions,
                 double excitement_factor);

#endif
