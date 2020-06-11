// Filthy awful hack to test.
#include "gtest/gtest.h"
#include <Eigen/Dense>
#include <memory>
#define private public
#include "velocity_verlet.h"
#include "energies.h"


TEST(VerletTest, StationaryEnergy)
{
    Eigen::VectorXd rail(2);
    rail << 1.0, 0.0;
    auto bond_ptr = std::make_unique<HarmonicBond>(1.0, 1.0, 1, 0, rail);
    std::vector<std::unique_ptr<Bond>> bonds_vec;
    bonds_vec.push_back(std::move(bond_ptr));
    
    Eigen::MatrixXd positions(2, 2);
    positions << 0.0, 0.0,
                 1.0, 0.0;
    auto forces = calculate_forces(bonds_vec, positions);
    for (int row = 0; row < forces.rows(); ++row)  {
        for (int col = 0; col < forces.cols(); ++col)  {
            ASSERT_FLOAT_EQ(forces(row, col), 0.0);
        }    
    }
}

TEST(VerletTest, ConservationOfEnergy)
{
    Eigen::VectorXd rail(2);
    rail << 1.0, 0.0;
    auto bond_ptr = std::make_unique<HarmonicBond>(1.0, 1.0, 1, 0, rail);
    std::vector<std::unique_ptr<Bond>> bonds_vec;
    bonds_vec.push_back(std::move(bond_ptr));
    
    Eigen::MatrixXd positions(2, 2);
    positions << 0.0, 0.0,
                 1.0, 0.0;
                 
    Eigen::MatrixXd velocities(2, 2);
    positions << 0.0, 0.0,
                 0.0, 0.0;
                 
    Eigen::MatrixXd accelerations(2, 2);
    positions << 0.0, 0.0,
                 0.0, 0.0;
                 
    Eigen::VectorXd masses(2);
    masses << 1.0, 1.0;
    
    excite_bond(bonds_vec[0], positions, masses, 1.5);             
    const auto initial_energy = calculate_kinetic_energy(velocities, masses) + calculate_potential_energy(bonds_vec, positions);
    
    constexpr int num_steps = 10000;
    constexpr double timestep = 0.01;
    constexpr double error_threshold = 1e-10;
    for (int i = 0; i < num_steps; ++i) {
        velocity_verlet(positions, velocities,
                        accelerations, masses,
                        timestep, bonds_vec);
        const auto total_energy = calculate_kinetic_energy(velocities, masses) + calculate_potential_energy(bonds_vec, positions);
        ASSERT_NEAR(initial_energy, total_energy, error_threshold);
    }
}
