// ============================================================================
// physical_constants.cpp - 物理常数库实现
// ============================================================================

#include "units/physical_constants.h"
#include <map>

namespace units {

const std::vector<PhysicalConstant>& PhysicalConstants::get_all() {
    static const std::vector<PhysicalConstant> kConstants = {
        {"c_light", "c", Quantity(299792458.0L, "m/s"), "Speed of light in vacuum"},
        {"h_planck", "h", Quantity(6.62607015e-34L, "J*s"), "Planck constant"},
        {"hbar", "hbar", Quantity(1.054571817e-34L, "J*s"), "Reduced Planck constant (h / 2pi)"},
        {"G_gravity", "G", Quantity(6.67430e-11L, "m^3/(kg*s^2)"), "Newtonian gravitational constant"},
        {"g_earth", "g0", Quantity(9.80665L, "m/s^2"), "Standard acceleration of gravity on Earth"},
        {"e_charge", "e_c", Quantity(1.602176634e-19L, "C"), "Elementary charge"},
        {"m_electron", "m_e", Quantity(9.1093837015e-31L, "kg"), "Electron rest mass"},
        {"m_proton", "m_p", Quantity(1.67262192369e-27L, "kg"), "Proton rest mass"},
        {"m_neutron", "m_n", Quantity(1.67492749804e-27L, "kg"), "Neutron rest mass"},
        {"k_boltzmann", "k_B", Quantity(1.380649e-23L, "J/K"), "Boltzmann constant"},
        {"N_avogadro", "N_A", Quantity(6.02214076e23L, "1/mol"), "Avogadro constant"},
        {"R_gas", "R", Quantity(8.314462618L, "J/(mol*K)"), "Molar gas constant"},
        {"eps_0", "eps_0", Quantity(8.8541878128e-12L, "F/m"), "Vacuum electric permittivity"},
        {"mu_0", "mu_0", Quantity(1.25663706212e-6L, "N/A^2"), "Vacuum magnetic permeability"},
        {"sigma_sb", "sigma_sb", Quantity(5.670374419e-8L, "W/(m^2*K^4)"), "Stefan-Boltzmann constant"},
        {"atm_pressure", "atm_p", Quantity(101325.0L, "Pa"), "Standard atmospheric pressure"}
    };
    return kConstants;
}

const PhysicalConstant* PhysicalConstants::find(const std::string& name_or_symbol) {
    const auto& all = get_all();
    for (const auto& c : all) {
        if (c.name == name_or_symbol || c.symbol == name_or_symbol) {
            return &c;
        }
    }
    return nullptr;
}

} // namespace units
