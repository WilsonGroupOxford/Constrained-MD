#include "velocity_verlet.h"

#include <Eigen/Dense>
#include <memory>
#include <vector>

#include "bond.h"
#include "constants.h"


Eigen::MatrixXd calculate_forces(const std::vector<std::unique_ptr<Bond>>& bonds,
                                 const Eigen::MatrixXd& positions) {
    //! Calculate the forces on each atom
    //! uses Newton's third law in that each force has an equal and opposite reaction
    /*! 
     *! \param[in] bonds a vector of polymorphic bonds between two atoms which will be used to calculate forces. Must have a force virtual method
     *! \param[in] positions the positions, used by the bonds to calculate forces
     */
    Eigen::MatrixXd forces = Eigen::MatrixXd::Zero(positions.rows(), positions.cols());
    for (const auto& bond : bonds) {
        Eigen::VectorXd bond_force = bond->force(positions);
        forces.row(bond->atoms[0]) -= bond_force;
        forces.row(bond->atoms[1]) += bond_force;
    }
    return forces;
}

Eigen::MatrixXd velocity_verlet(Eigen::MatrixXd& positions, Eigen::MatrixXd& velocities,
                                Eigen::MatrixXd& accelerations, const Eigen::VectorXd& masses,
                                const double dt, const std::vector<std::unique_ptr<Bond>>& bonds) {
    //! Solve Newton's laws of motion numerically by stepping forwards a small timestep.
    //! Solve $ F = m a $ numerically in increments of $ \delta t $. Use the equations
    //! $ x = x + v * \delta t + 0.5 * a * \delta t^2 $
    //! $ a_new = F / m $
    //! $ v = v + 0.5 * (a_new + a_old) $
    /*! 
     *! \param[in] positions an array of particle positions which we will mutate
     *! \param[in] velocities an array of particle velocities which we will mutate
     *! \param[in] accelerations  an array of particle accelerations which we will mutate
     *! \param[in] masses an array of particle masses
     *! \param[in] dt the timestep, which should be at least 1/10 of the period of the shortest bond frequency
     *! \param[in] bonds an array of bonds which will be used to calculate forces
     *! \return new_positions the new updated positions
     */
    positions += velocities * dt + (0.5 * accelerations * dt * dt);
    Eigen::MatrixXd new_accelerations = calculate_forces(bonds, positions);
    // TODO: Speed this up by using Eigen array
    for (int i = 0; i < new_accelerations.rows(); ++i) {
        new_accelerations.row(i) /= masses[i];
    }
    velocities += 0.5 * (accelerations + new_accelerations) * dt;
    accelerations = new_accelerations;
    return positions;
}

void excite_bond(const std::unique_ptr<Bond>& bond, Eigen::MatrixXd& positions,
                 double excitement_factor) {
    //! Excite a single bond by stretching it by the excitement factor.
    //! For the bond, replace the position of the atoms to be exactly along the bond
    //! vector, with the distance between them being (equilibrium distance * excitement factor).
    /*! 
     *! \param[in] bond a bond object to excite
     *! \param[in] positions a mutable array of positions, which we will change
     *! \param[in] excitement_factor the multiple of the equilibrium bond length to set the new bond length to
     */
    const Eigen::VectorXd displacement_vector
        = bond->rail * bond->equilibrium_distance * excitement_factor * -1;
    const Eigen::VectorXd old_position = positions.row(bond->atoms[0]);
    positions.row(bond->atoms[1]) = old_position + displacement_vector;
}
