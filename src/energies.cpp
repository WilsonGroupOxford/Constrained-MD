#include "energies.h"

#include "bond.h"

#include <Eigen/Dense>
#include <memory>
#include <vector>
#include <algorithm>
#include <numeric>


double calculate_kinetic_energy(const Eigen::MatrixXd& velocities, const Eigen::MatrixXd& masses) {
    return (masses.transpose()  * velocities.rowwise().squaredNorm()).sum() / 2.0;
}

double calculate_potential_energy(const std::vector<std::unique_ptr<Bond>>& bonds, const Eigen::MatrixXd& positions) {
    std::vector<double> potential_energies;
    potential_energies.reserve(bonds.size());
    std::transform(bonds.begin(), bonds.end(), std::back_inserter(potential_energies),
    [positions](const auto& bond){ return bond->energy(positions);});
    
    return std::accumulate(potential_energies.begin(),
                                                  potential_energies.end(),
                                                    0.0
                                                   );
}
