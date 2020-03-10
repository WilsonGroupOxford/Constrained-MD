#include <Eigen/Dense>
#include <fstream>
#include <iostream>
#include <vector>

class Bond {
public:
    const double force_constant;
    const double equilibrium_distance;
    const int first_atom;
    const int second_atom;
    const Eigen::Vector2d rail;
    Bond(double force_constant, double equilibrium_distance, int first_atom, int second_atom,
         Eigen::Vector2d rail)
        : force_constant(force_constant)
        , equilibrium_distance(equilibrium_distance)
        , first_atom(first_atom)
        , second_atom(second_atom)
        , rail(rail.normalized()) {};

    Eigen::Vector2d project_onto_rail(const Eigen::Vector2d& input) const {
        return input.dot(rail) * rail;
    }
    Eigen::Vector2d force(const Eigen::MatrixXd& positions) const {
        Eigen::Vector2d separation = positions.row(first_atom) - positions.row(second_atom);
        Eigen::Vector2d force_direction = force_constant * (separation - (rail * equilibrium_distance));
        return project_onto_rail(force_direction);
    }

    double energy(const Eigen::MatrixXd& positions) const {

        double distance = (positions.row(first_atom) - positions.row(second_atom)).norm();
        double distance_from_eqm = (distance - equilibrium_distance);
        return distance_from_eqm * distance_from_eqm * force_constant;
    }
};

Eigen::MatrixXd normalise_vectors(Eigen::MatrixXd& vectors) {
    return vectors.rowwise().normalized();
}

Eigen::MatrixXd calculate_forces(const std::vector<Bond>& bonds, const Eigen::MatrixXd& positions) {
    Eigen::MatrixXd forces = Eigen::MatrixXd::Zero(positions.rows(), positions.cols());
    for (const auto& bond : bonds) {
        Eigen::Vector2d bond_force = bond.force(positions);
        forces.row(bond.first_atom) -= bond_force;
        forces.row(bond.second_atom) += bond_force;
    }
    return forces;
}

Eigen::MatrixXd velocity_verlet(Eigen::MatrixXd& positions, Eigen::MatrixXd& velocities,
                                Eigen::MatrixXd& accelerations, const double dt,
                                const std::vector<Bond>& bonds) {
    positions += velocities * dt + (0.5 * accelerations * dt * dt);
    Eigen::MatrixXd new_accelerations = calculate_forces(bonds, positions);
    velocities += 0.5 * (accelerations + new_accelerations) * dt;
    accelerations = new_accelerations;
    return positions;
}

void write_xyz(const std::string& filename, const Eigen::MatrixXd& positions,
               const Eigen::MatrixXd& velocities, const Eigen::MatrixXd& accelerations,
               const int step) {
    auto output_file = std::ofstream(filename, std::ios::app);

    output_file << positions.rows() << "\n";
    output_file << "# FRAME " << step << "\n";
    std::array ELEMENTS { "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne" };
    for (int i = 0; i < positions.rows(); ++i) {
        output_file << ELEMENTS[i] << "\t" << positions.row(i) << "\t0.0\t" << velocities.row(i)
                    << "\t0.0\t" << accelerations.row(i) << "\t 0.0 \t"
                    << "\n";
    }
}

Eigen::MatrixXd load_positions(const std::string& filename) {
    std::ifstream input_file {filename};
    if (! input_file.good()) {
        throw std::runtime_error("Could not open " + filename);
    }

    std::vector<Eigen::Vector2d> positions_vec;
    while (true) {
        std::string line;
        std::getline(input_file, line);
        if (input_file.eof()) {
            break;
        }
        std::cout << line << "\n";
        Eigen::Vector2d pos_vec;
        std::istringstream iss{line};
        iss >> pos_vec[0];
        iss >> pos_vec[1];
        if (iss.fail()) {
            throw std::runtime_error("Could not convert " + line);
        }
        positions_vec.push_back(std::move(pos_vec));
    }
    Eigen::MatrixXd positions(positions_vec.size(), 2);
    for (int i = 0; i < static_cast<int>(positions_vec.size()); ++i) {
        positions.row(i) = std::move(positions_vec[i]);
    }
    return positions;
}

std::vector<Bond> load_bonds(const std::string& filename, const Eigen::ArrayXXd& positions) {
    //! Load bond data from a file.

    std::vector<Bond> bonds;
    std::ifstream input_file{filename};
    if (! input_file.good()) {
        throw std::runtime_error("Could not open " + filename);
    }

    while (true) {
        std::string line;
        std::getline(input_file, line);
        if (input_file.eof()) {
            break;
        }
        double equilibrium_distance;
        double force_constant;
        int first_atom;
        int second_atom;

        std::istringstream iss{line};
        iss >> equilibrium_distance;
        iss >> force_constant;
        iss >> first_atom;
        iss >> second_atom;

        bonds.emplace_back(equilibrium_distance, force_constant, first_atom, second_atom,
                           positions.row(first_atom) - positions.row(second_atom));
    }
    return bonds;
}
int main() {
    {
        Eigen::MatrixXd positions = load_positions("positions.dat");
        auto bonds = load_bonds("bonds.dat", positions);

        Eigen::MatrixXd velocities = Eigen::MatrixXd::Zero(positions.rows(), positions.cols());
        Eigen::MatrixXd accelerations = Eigen::MatrixXd::Zero(positions.rows(), positions.cols());
        std::ofstream output_file { "test.xyz" };
        for (int step = 0; step < 10000; ++step) {
            velocity_verlet(positions, velocities, accelerations, 0.01, bonds);
            write_xyz("test.xyz", positions, velocities, accelerations, step);
        }
    }
    return 0;
}
