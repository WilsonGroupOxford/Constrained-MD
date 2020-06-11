#include <Eigen/Dense>
#include <fstream>
#include <iostream>
#include <memory>
#include <vector>
#include <string_view>
#include <numeric>

#include "yaml-cpp/yaml.h"

#include "bond.h"
#include "constants.h"
#include "file_io.h"
#include "velocity_verlet.h"
#include "energies.h"

const std::string DEFAULT_CONFIG_FILE{"config.yaml"};

template <typename BondIter>
bool check_bond_sampling(BondIter bonds_begin, BondIter bonds_end, const Eigen::MatrixXd& masses, const double timestep, const float min_samples=10) {
    std::vector<double> bond_periods;
    std::transform(bonds_begin, bonds_end, std::back_inserter(bond_periods), [masses](auto&& bond){return bond->period(masses);});
    const double min_bond_period = *std::min_element(bond_periods.begin(), bond_periods.end());
    return (timestep * min_samples < min_bond_period);
}



int main(int argc, char** argv) {
    std::string config_filename;
    if (argc > 1) {
        config_filename = argv[1];
    } else {
        config_filename = DEFAULT_CONFIG_FILE;
    }
    YAML::Node config;
    try {
         config = YAML::LoadFile(config_filename);
    } catch (const YAML::BadFile& ex) {
         std::cerr << "Could not open " << config_filename << ". Is it present in this directory?\n";
         throw ex;
    }
    const std::string positionfile = config["positionfile"].as<std::string>();
    const std::string bondfile = config["bondfile"].as<std::string>();
    const std::string outputfile = config["outputfile"]["filename"].as<std::string>();
    const int output_frequency = config["outputfile"]["frequency"].as<int>();
    const int console_frequency = config["outputfile"]["consolefrequency"].as<int>();
    const double timestep = config["simulation"]["timestep"].as<double>();
    const int number_steps = config["simulation"]["steps"].as<int>();

    constexpr auto unit_type = UnitType::ARBITRARY;
    auto [positions, atom_names] = load_positions(positionfile, unit_type, 3);
    auto bonds = load_bonds(bondfile, positions, atom_names);

    Eigen::VectorXd masses = load_masses(atom_names, unit_type);
    if (!check_bond_sampling(bonds.begin(), bonds.end(), masses, timestep, 10)) {
        std::cerr << "Shortest bond period is sampled less than ten times per timestep.\n";
    }

    Eigen::MatrixXd velocities = Eigen::MatrixXd::Zero(positions.rows(), positions.cols());
    

    for (const auto& excitement : config["excitements"]) {
        const auto atom_a_name = excitement["atoms"][0].as<std::string>();
        const int atom_a_id = atom_name_to_id(atom_names.begin(), atom_names.end(), atom_a_name);

        const auto atom_b_name = excitement["atoms"][1].as<std::string>();
        const int atom_b_id = atom_name_to_id(atom_names.begin(), atom_names.end(), atom_b_name);

        // Now find the correct bond.
        auto correct_bond = std::find_if(bonds.begin(), bonds.end(), [atom_a_id, atom_b_id](const auto& bond){
            return (bond->atoms[0] == atom_a_id) && (bond->atoms[1] == atom_b_id);
        });
        excite_bond(*correct_bond, positions, masses, excitement["factor"].as<double>());

        std::cout << "Exciting a bond between " << atom_a_name << " and " << atom_b_name << "\n";
    }
    
    // We need to set the initial accelerations to be correct.
    Eigen::MatrixXd accelerations = calculate_accelerations(bonds, positions, masses);

    std::ofstream output_file { outputfile };
    std::ofstream bond_excitement_file { "bonds.csv" };
    
    bond_excitement_file << "Time, ";
    for (const auto& bond : bonds) {
        bond_excitement_file << atom_names[bond->atoms[0]] << "->" << atom_names[bond->atoms[1]]
                             << ",";
    }
    bond_excitement_file << "\n";
    for (int step = 0; step < number_steps; ++step) {
        velocity_verlet(positions, velocities, accelerations, masses, timestep, bonds);
        if (step % console_frequency == 0 ) {;
            const auto potential_energy = calculate_potential_energy(bonds, positions);
            const auto kinetic_energy = calculate_kinetic_energy(velocities, masses);
            std::cout << "Step " << step << ", time = " << step * timestep << ", PE = " << potential_energy << ", KE = "
<< kinetic_energy << " TE = " << potential_energy + kinetic_energy  << "\n";
        }
        if (step % output_frequency == 0) {
            write_xyz(outputfile, positions, velocities, accelerations, step, atom_names);
            bond_excitement_file << step * timestep << ", ";
            for (const auto& bond : bonds) {
                bond_excitement_file << bond->get_excitement_factor(positions)<< ",";
            }
            bond_excitement_file << "\n";
        }
    }
    return 0;
}
