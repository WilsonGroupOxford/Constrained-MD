#ifndef BOND_H
#define BOND_H

#include <Eigen/Dense>
#include <array>
#include <algorithm>
#include <string_view>
#include <iterator>
#include <sstream>

template <typename Iter>
constexpr int atom_name_to_id(Iter atom_names_begin, Iter atom_names_end, std::string_view atom_name) {
    //! Convert an atom name to a numerical index.
    //! Searches through the list of names to find an exactly equal name,
    //! and return the distance from atom_names_begin.  Returns the first index
    //! if there are duplicate names (which there should not be!), and throws
    //! a runtime error if we cannot find the name.
    /*!
     *! \param[in] atom_names_begin the iterator to the start of the atom names array
     *! \param[in] atom_names_end the iterator to one-past-the-end of the atom names array
     *! \param[in] atom_name the name of the atom to find
     *! \return index the index offset of atom_name in the array from atom_names_begin. If atom_names_begin is atom_names.begin(), this is the index of the atom to access.
     */
    const auto name_pos = std::find(atom_names_begin, atom_names_end, atom_name);
    if (name_pos == atom_names_end) {
        const char* const delim = ", ";

        std::ostringstream imploded;
        std::copy(atom_names_begin, atom_names_end,
                   std::ostream_iterator<std::string>(imploded, delim));

        throw std::runtime_error("Could not find " + std::string(atom_name) + " in the list of atom names." +
                                 "Valid names are: " + imploded.str());
    }
    return std::distance(
            atom_names_begin, name_pos);
}

enum class BondType {
    HARMONIC,
};

class Bond {
public:
    const double force_constant;
    const double equilibrium_distance;
    const std::array<int, 2> atoms;
    const Eigen::VectorXd rail;
    Bond(double force_constant, double equilibrium_distance, int first_atom, int second_atom,
         Eigen::VectorXd rail);

    Eigen::VectorXd project_onto_rail(const Eigen::VectorXd& input) const;
    double get_excitement_factor(const Eigen::MatrixXd& positions) const;

    virtual Eigen::VectorXd force(const Eigen::MatrixXd& positions) const = 0;
    virtual double energy(const Eigen::MatrixXd& positions) const = 0;
    virtual double frequency(const Eigen::VectorXd& masses) const = 0;
    virtual double period(const Eigen::VectorXd& masses) const = 0;
};

class HarmonicBond : public Bond {
public:
    HarmonicBond(double force_constant, double equilibrium_distance, int first_atom,
                 int second_atom, Eigen::VectorXd rail);
    Eigen::VectorXd force(const Eigen::MatrixXd& positions) const override;

    double energy(const Eigen::MatrixXd& positions) const override;
    double frequency(const Eigen::VectorXd& masses) const override;
    double period(const Eigen::VectorXd& masses) const override;
};
#endif // BOND_H
