#include "file_io.h"

#include "bond.h"
#include "constants.h"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <memory>
#include <set>
#include <string>
#include <tuple>
#include <vector>
namespace {

bool str_starts_with(const std::string& input, const std::string& prefix){
    //! Overload for missing standard library functions while we wait for C++20
    return input.rfind(prefix, 0) == 0;
}

template <typename Container>
constexpr bool str_contains(const Container& input, const std::string& substr){
    //! Overload for missing standard library functions while we wait for C++20
    return input.find(substr) != input.end();
}

std::unique_ptr<HarmonicBond> construct_harmonic_bond(std::istringstream& iss,
                                                      const Eigen::ArrayXXd& positions,
                                                      const std::vector<std::string>& atom_names,
                                                      const UnitType& unit_type) {
    //! Construct  a Harmonic Bond from a string.
    //! Uses a stringstream representing one line from a file to construct a bond.
    //! The stringstream should be in the following format:
    //! first_atom second_atom force_constant equilibrium_distance
    //! This looks up the atom names in the atom_names vector to get their IDs
    //! for efficient lookups in the bonds.
    /*! \param[in] iss a stringstream containing one line of a file
     *! \param[in] positions the positions of the atoms to construct rails from
     *! \param[in] atom_names an vector of atom names that we can look them up in
     *! \param[in] unit_type the units we will use for this computation, and we adjust those into atomic units
     *! \return a pointer to a Harmonic Bond disguised as a base class Bond
     */
    std::cout << iss.str() << "\n";
    double equilibrium_distance;
    double force_constant;
    std::string first_atom;
    std::string second_atom;

    iss >> first_atom;
    iss >> second_atom;
    iss >> force_constant;
    iss >> equilibrium_distance;

    int first_atom_id = atom_name_to_id(atom_names.begin(), atom_names.end(), first_atom);
    int second_atom_id = atom_name_to_id(atom_names.begin(), atom_names.end(), second_atom);

    switch (unit_type) {
    case UnitType::ARBITRARY:
        break;
    case UnitType::ATOMIC:
        // Convert the distance from angstroms into bohr radii,
        // and the force constant from wavenumbers into hartrees
        equilibrium_distance /= constants::bohr_ang;
        force_constant = force_constant * constants::hartree / constants::wavenumber;
        break;
    case UnitType::REAL:
        equilibrium_distance /= constants::bohr;
        force_constant /= constants::hartree;
        break;
    }
    Eigen::VectorXd separation = positions.row(first_atom_id) - positions.row(second_atom_id);
    return std::make_unique<HarmonicBond>(std::move(force_constant),
                                          std::move(equilibrium_distance), std::move(first_atom_id),
                                          std::move(second_atom_id), std::move(separation));
}

std::string remove_numbers(const std::string& str) {
    //! Remove numbers from a string.
    /*!
     * \param str the string to remove numbers from
     * \return a new string with all numbers removed
     */
    auto str_copy = str;
    str_copy.erase(
                remove_if(str_copy.begin(), str_copy.end(), [](char c) { return std::isdigit(c); }),
                str_copy.end());
    return str_copy;
}

std::string atom_name_to_element(const std::string& input, const char comment = '%',
                                 const std::set<std::string>& valid_elements
                                 = constants::element_symbols) {
    //! Convert an atom name to an element
    //! Split the atom name at the comment symbol, by default "%",
    //! and strip out all numbers. Return the element symbol.
    /*!
         * \param input an input in the form "Ab%comment123".
         * \return element a one or-two letter code verified as a correct element
         */
    std::string elem = input.substr(0, input.find(comment, 0));
    if (elem == input) {
        elem = remove_numbers(input);
    }

    if (!str_contains(valid_elements, elem)) {
        throw std::runtime_error("Could not convert " + input + " to a valid element.");
    }
    return elem;
}

}

std::tuple<Eigen::MatrixXd, std::vector<std::string>> load_positions(const std::string& filename,
                                                                     const UnitType& unit_type,
                                                                     const int dimension
                                                                     ) {
    /*!
     * \brief Load positions from an .xyz file
     //! Files lines in the format
     //! atom_name position_x position_y position_z
     //! First line is the number of atoms, second is the comment which is ignored.
     *! \param filename the name of the file to read in from
     *! \param unit_type the units to adjust these into
     *! \param dimension can be 2D or 3D, represents the dimensionality of the simulation
     */
    if (dimension != 2 && dimension != 3) {
        throw std::runtime_error("Dimension must be 2 or 3. Instead we got " + std::to_string(dimension));
    }
    std::ifstream input_file { filename };
    if (!input_file.good()) {
        throw std::runtime_error("Could not open " + filename);
    }

    std::vector<Eigen::VectorXd> positions_vec;
    std::vector<std::string> atom_types;
    std::string dummy;
    // Number of atoms

    std::getline(input_file, dummy);
    int num_atoms = std::stoi(dummy);
    // Comment
    std::getline(input_file, dummy);
    while (true) {
        std::string line;
        std::getline(input_file, line);
        if (input_file.eof()) {
            break;
        }
        Eigen::VectorXd pos_vec = Eigen::VectorXd::Zero(dimension);
        std::istringstream iss { line };
        std::string atom_type;
        iss >> atom_type;
        atom_types.push_back(atom_type);
        for (int i = 0; i < dimension; ++i) {
            iss >> pos_vec[i];
            std::cout << i << ", " <<pos_vec[i] << ", ";
        }
        std::cout << "\n";
        if (iss.fail()) {
            throw std::runtime_error("Could not convert " + line);
        }
        switch (unit_type) {
        case UnitType::ARBITRARY:
            // Don't do anything to the vectors
            break;
        case UnitType::ATOMIC:
            pos_vec /= constants::bohr_ang;
            break;
        case UnitType::REAL:
            pos_vec /= constants::bohr;
            break;
        }
        positions_vec.push_back(std::move(pos_vec));
    }
    if (static_cast<int>(positions_vec.size()) != num_atoms) {
        throw std::runtime_error("Number of positions read (" + std::to_string(positions_vec.size()) +
                                 ") does not equal declared number of atoms (" + std::to_string(num_atoms) + ")");
    }
    Eigen::MatrixXd positions(positions_vec.size(), dimension);
    for (int i = 0; i < static_cast<int>(positions_vec.size()); ++i) {
        positions.row(i) = std::move(positions_vec[i]);
    }

    auto unique_it = std::unique(atom_types.begin(), atom_types.end());
    if (unique_it != atom_types.end()) {
        throw std::runtime_error("Warning: Atom names are not unique.");
    }
    return { positions, atom_types };
}

std::vector<std::unique_ptr<Bond>> load_bonds(const std::string& filename,
                                              const Eigen::ArrayXXd& positions,
                                              const std::vector<std::string>& atom_names) {
    //! Load bond data from a file.
    //!
    //! Read in from filename a list of harmonic bonds, with equilibrium distances,
    //! force constants and two atoms that they are between. Each bond can only
    //! apply force along its original vector, so calculate those as well
    //! from positions. Convert the bond lengths and distances into atomic units.
    //! TODO 2020-03-11(MB) -- Convert this to a better file format?
    //!
    /*!
     *! \param filename the name of the file to read from..
     *! \param positions an array of atomic positions to calculate bond rails from.
     *! \param unit_type an enum dictating the unit system of the file.
     *! \return bonds an iterable of all the bonds in the system.
     */

    int lines_read = 0;
    std::vector<std::unique_ptr<Bond>> bonds;
    bonds.reserve(positions.rows());
    std::ifstream input_file { filename };
    if (!input_file.good()) {
        throw std::runtime_error("Could not open " + filename);
    }
    std::string units_line;
    std::getline(input_file, units_line);
    ++lines_read;
    UnitType unit_type;
    if (str_starts_with(units_line, "atomic")) {
        unit_type = UnitType::ATOMIC;
    } else if (str_starts_with(units_line, "real")) {
        unit_type = UnitType::REAL;
    } else if (str_starts_with(units_line, "arbitrary"))  {
        unit_type = UnitType::ARBITRARY;
    } else {
        throw std::runtime_error("Line " + std::to_string(lines_read) + ":Could not convert " + units_line + " to a units type. Expected: atomic, real, arbitrary");
    }
   

    while (true) {
        std::string line;
        std::getline(input_file, line);
        ++lines_read;
        if (input_file.eof()) {
            break;
        }

        std::istringstream iss { line };
        std::string bond_type_str;
        iss >> bond_type_str;
        BondType bond_type;
        if (bond_type_str == "harmonic") {
            bond_type = BondType::HARMONIC;
        } else {
            throw std::runtime_error("Line " + std::to_string(lines_read) + ": Could not convert " + bond_type_str + " to a bond type.");
        }

        switch (bond_type) {
        case BondType::HARMONIC:
            std::unique_ptr<Bond> this_bond
                    = construct_harmonic_bond(iss, positions, atom_names, unit_type);
            bonds.push_back(std::move(this_bond));
            break;
        }
    }
    return bonds;
}

Eigen::VectorXd load_masses(const std::vector<std::string>& atom_types, const UnitType& unit_type) {
    //! Calculate masses from the types of atoms in the system.
    //! Use the elemental masses from the atom names.
    /*!
     *! \param atom_types the names of each atom in the system
     *! \param unit_type the unit type we use, does not convert if "arbitrary" and converts for others
     *! \return masses an eigen array of masses in the correct format
     */
    std::vector<double> atom_masses;
    atom_masses.reserve(atom_types.size());
    std::transform(atom_types.begin(), atom_types.end(), std::back_inserter(atom_masses),
                   [](auto atom_name) {
        return constants::element_masses.at(atom_name_to_element(atom_name));
    });
    Eigen::VectorXd eigen_atom_masses
            = Eigen::VectorXd::Map(atom_masses.data(), atom_masses.size());
    switch (unit_type) {
    case UnitType::ARBITRARY:
        break;
    case UnitType::ATOMIC:
        eigen_atom_masses = eigen_atom_masses * constants::dalton / constants::e_mass;
        break;
    case UnitType::REAL:
        eigen_atom_masses /= constants::e_mass;
        break;
    }

    return eigen_atom_masses;
}

void write_xyz(const std::string& filename, const Eigen::MatrixXd& positions,
               const Eigen::MatrixXd& velocities, const Eigen::MatrixXd& accelerations,
               const int step, const std::vector<std::string>& atom_types) {
    /*!
     * \brief Output the positions, accelerations and velocities to an xyz file.
     //! Appends to an existing file, so be sure to blank it out at the start
     //! of the run. Files lines in the format
     //! atom_name position_x position_y position_z velocity_x v_y v_z a_x a_y a_z
     //! Comment line is just the frame number
     *! \param filename the name of the file to output to
     *! \param positions an eigen array of atomic positions
     *! \param velocities an eigen array of atomic velocities
     *! \param accelerations an eigen array of accelerations
     *! \param step the step of the simulation this is
     *! \param atom_types the names to prefix the atoms with in the file
     */
    auto output_file = std::ofstream(filename, std::ios::app);

    output_file << positions.rows() << "\n";
    output_file << "# FRAME " << step << "\n";
    for (int i = 0; i < positions.rows(); ++i) {
        output_file << atom_types[i] << "\t" << positions.row(i) << "\t0.0\t" << velocities.row(i)
                    << "\t0.0\t" << accelerations.row(i) << "\t 0.0 \t"
                    << "\n";
    }
}
