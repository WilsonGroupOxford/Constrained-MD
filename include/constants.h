#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <algorithm>
#include <set>
#include <unordered_map>
enum class UnitType { ARBITRARY, ATOMIC, REAL };

namespace constants {
// Fundamental units
constexpr double e_mass = 9.1093837015e-31; // kilograms
constexpr double e_charge = 1.602176634e-19; // Coulomb
constexpr double planck = 6.62607015e-34; // Joule-second
constexpr double h_bar = 1.054571817e-34; // Joule-second
constexpr double inv_coulomb_constant = 1.11265005545e-10; // Farad per metre
constexpr double light_speed = 299792458; // metre per second
constexpr double boltz = 1.380649e-23; // Joule per kelvin
constexpr double avogadro = 6.02214076e23; // mol^-1

// Calculated units
constexpr double fine_structure
    = e_charge * e_charge / (inv_coulomb_constant * h_bar * light_speed);
constexpr double hartree
    = fine_structure * fine_structure * e_mass * light_speed * light_speed; // Joule
constexpr double bohr = h_bar / (e_mass * light_speed * fine_structure); // Metre
constexpr double bohr_ang = 1e10 * h_bar / (e_mass * light_speed * fine_structure);

constexpr double wavenumber = 1e-2 / (planck * light_speed); // cm^-1 -> Joule
// Conversion units out of atomic units

constexpr double dalton = 1e-3 / avogadro; // kg
namespace atomic {
    constexpr double time = h_bar / hartree; // seconds
    constexpr double velocity = bohr * hartree / h_bar; // Metres per second
    constexpr double force = hartree / bohr; // Newton
    constexpr double pressure = hartree / (bohr * bohr * bohr);
}

const std::unordered_map<std::string, double> element_masses = {
    { "H", 1.008 },        { "He", 4.002602 },     { "Li", 6.94 },     { "Be", 9.0121831 },
    { "B", 10.81 },        { "C", 12.011 },        { "N", 14.007 },    { "O", 15.999 },
    { "F", 18.998403163 }, { "Ne", 20.1797 },      { "Mg", 24.305 },   { "Al", 26.9815384 },
    { "Si", 28.085 },      { "P", 30.973761998 },  { "S", 32.06 },     { "Cl", 35.45 },
    { "Ar", 39.948 },      { "K", 39.0983 },       { "Ca", 40.078 },   { "Sc", 44.955908 },
    { "Ti", 47.867 },      { "V", 50.9415 },       { "Cr", 51.9961 },  { "Mn", 54.938043 },
    { "Fe", 55.845 },      { "Co", 58.933194 },    { "Ni", 58.6934 },  { "Cu", 63.546 },
    { "Zn", 65.38 },       { "Ga", 69.723 },       { "Ge", 72.63 },    { "As", 74.921595 },
    { "Se", 78.971 },      { "Br", 79.904 },       { "Kr", 83.798 },   { "Rb", 85.4678 },
    { "Sr", 87.62 },       { "Y", 88.90584 },      { "Zr", 91.224 },   { "Nb", 92.90637 },
    { "Mo", 95.95 },       { "Tc", 97 },           { "Ru", 101.07 },   { "Rh", 102.90549 },
    { "Pd", 106.42 },      { "Ag", 107.8682 },     { "Cd", 112.414 },  { "In", 114.818 },
    { "Sn", 118.71 },      { "Sb", 121.76 },       { "Te", 127.6 },    { "I", 126.90447 },
    { "Xe", 131.293 },     { "Cs", 132.90545196 }, { "Ba", 137.327 },  { "La", 138.90547 },
    { "Ce", 140.116 },     { "Pr", 140.90766 },    { "Nd", 144.242 },  { "Pm", 145 },
    { "Sm", 150.36 },      { "Eu", 151.964 },      { "Gd", 157.25 },   { "Tb", 158.925354 },
    { "Dy", 162.5 },       { "Ho", 164.930328 },   { "Er", 167.259 },  { "Tm", 168.934218 },
    { "Yb", 173.045 },     { "Lu", 174.9668 },     { "Hf", 178.486 },  { "Ta", 180.94788 },
    { "W", 183.84 },       { "Re", 186.207 },      { "Os", 190.23 },   { "Ir", 192.217 },
    { "Pt", 195.084 },     { "Au", 196.96657 },    { "Hg", 200.592 },  { "Tl", 204.38 },
    { "Pb", 207.2 },       { "Bi", 208.9804 },     { "Po", 209 },      { "At", 210 },
    { "Rn", 222 },         { "Fr", 223 },          { "Ra", 226 },      { "Ac", 227 },
    { "Th", 232.0377 },    { "Pa", 231.03588 },    { "U", 238.02891 }, { "Np", 237 },
    { "Pu", 244 },         { "Am", 243 },          { "Cm", 247 },      { "Bk", 247 },
    { "Cf", 251 },         { "Es", 252 },          { "Fm", 257 },      { "Md", 258 },
    { "No", 259 },         { "Lr", 262 },          { "Rf", 267 },      { "Db", 270 },
    { "Sg", 269 },         { "Bh", 270 },          { "Hs", 270 },      { "Mt", 278 },
    { "Ds", 281 },         { "Rg", 281 },          { "Cn", 285 },      { "Nh", 286 },
    { "Fl", 289 },         { "Mc", 289 },          { "Lv", 293 },      { "Ts", 293 },
    { "Og", 294 },

};

template <typename T, typename U>
constexpr std::set<T> get_unordered_map_keys(const std::unordered_map<T, U> & map) {
    std::set<T> new_set;
    std::transform(map.begin(), map.end(),
                   std::inserter(new_set, new_set.begin()),
                   [](auto pair) { return pair.first; });
    return new_set;
}

const auto element_symbols = get_unordered_map_keys(element_masses);
}



#endif // CONSTANTS_H
