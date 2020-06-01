#include "bond.h"
#include <Eigen/Dense>
#include <iostream>

Bond::Bond(double force_constant, double equilibrium_distance, int first_atom, int second_atom,
           Eigen::VectorXd rail)
    : force_constant(force_constant)
    , equilibrium_distance(equilibrium_distance)
    , atoms { first_atom, second_atom }
    , rail {rail.normalized()} {}

Eigen::VectorXd Bond::project_onto_rail(const Eigen::VectorXd& input) const {
    //! Project a given vector onto this bond rail
    /*!
     *! \param input the vector to do the projection of
     *! \return projected the component of the vector in the direction of the rail
     */
    return input.dot(rail) * rail;
}

double Bond::get_excitement_factor(const Eigen::MatrixXd& positions) const {
    Eigen::VectorXd separation = positions.row(atoms[0]) - positions.row(atoms[1]);
    separation = project_onto_rail(separation);
    double distance = separation.dot(rail);
    return distance / equilibrium_distance;

}

HarmonicBond::HarmonicBond(double force_constant, double equilibrium_distance, int first_atom,
                           int second_atom, Eigen::VectorXd rail)
    : Bond(force_constant, equilibrium_distance, first_atom, second_atom, rail) {}

Eigen::VectorXd HarmonicBond::force(const Eigen::MatrixXd& positions) const {
    const Eigen::VectorXd separation = positions.row(atoms[0]) - positions.row(atoms[1]);
    const Eigen::VectorXd force_direction
            = force_constant * (separation - (rail * equilibrium_distance));
    return project_onto_rail(force_direction);
}

double HarmonicBond::frequency(const Eigen::VectorXd& masses) const {
    const auto reduced_mass = masses[atoms[0]] * masses[atoms[1]] / (masses[atoms[0]] + masses[atoms[1]]);
    return std::sqrt(force_constant / reduced_mass);
}

double HarmonicBond::period(const Eigen::VectorXd& masses) const {
    return 2 * M_PI / frequency(masses);
}

double HarmonicBond::energy(const Eigen::MatrixXd& positions) const {
    double distance = project_onto_rail(positions.row(atoms[0]) - positions.row(atoms[1])).norm();
    double distance_from_eqm = (distance - equilibrium_distance);
    return distance_from_eqm * distance_from_eqm * force_constant;
}
