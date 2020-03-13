#ifndef BOND_H
#define BOND_H

#include <Eigen/Dense>
#include <array>

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
};

class HarmonicBond : public Bond {
public:
    HarmonicBond(double force_constant, double equilibrium_distance, int first_atom,
                 int second_atom, Eigen::VectorXd rail);
    Eigen::VectorXd force(const Eigen::MatrixXd& positions) const override;

    double energy(const Eigen::MatrixXd& positions) const override;
};
#endif // BOND_H
