#ifndef _CONSTANTS_H
#define _CONSTANTS_H
#include <cmath>

namespace celerium{
namespace phys{

const double fine_structure = 0.0072973525664;
const double one_over_fine_structure = 137.035999139;

const double planck_SI  = 6.626070040e-34; 
const double planck_AeV = 4.135667662e-15;
const double planck_He  = 2.0*M_PI; 
const double planck_Ry  = 2.0*M_PI;

const double planck_bar_SI  = 1.054571800e-34;
const double planck_bar_AeV = 6.582119514e-16;
const double planck_bar_He  = 1.0; 
const double planck_bar_Ry  = 1.0;

const double light_speed_SI = 2.99792458e8;
const double light_speed_AeV = 2.99792458e18;
const double light_speed_He = one_over_fine_structure;
const double light_speed_Ry = 2.0*one_over_fine_structure;

// C = K e^2 from Coulomb's law, C = alpha hbar c
const double coulomb_coeff_SI  = 2.3070775128274957426e-28;
const double coulomb_coeff_AeV = 14.3996453513098353096;
const double coulomb_coeff_He  = 1.0;
const double coulomb_coeff_Ry  = 2.0;

} // end of namespace physics

namespace convert{

// SI to AeV:
const double  m_to_A = 1e10;
const double J_to_eV = 1.0/1.6021766208e-19;

// SI to He:
const double m_to_a0 = 1.0/5.2917721067e-11;
const double J_to_He = 1.0/4.35974417e-18; 

// SI to Ry:
//const double m_to_a0 = 1.0/5.2917721067e-11;
const double J_to_Ry = 2.0/4.35974417e-18;

// AeV to SI:
const double  A_to_m = 1e-10;
const double eV_to_J = 1.6021766208e-19;

// AeV to He:
const double  A_to_a0 = 1.0/5.2917721092e-1;
const double eV_to_He = 0.5/13.605693009;

// AeV to Ry:
//const double  A_to_a0 = 1.0/5.2917721092e-1;
const double eV_to_Ry = 1.0/13.605693009; 

// He to SI:
const double a0_to_m = 5.2917721092e-11;
const double He_to_J = 4.35974417e-18;

// He to AeV:
const double  a0_to_A  = 5.2917721092e-1;
const double He_to_eV = 2.0*13.605693009;

// He to Ry
const double He_to_Ry = 2.0;

// Ry to SI:
//const double a0_to_m = 5.2917721092e-11;
const double Ry_to_J = 0.5*4.35974417e-18;

// Ry to AeV:
//const double  a0_to_A  = 5.2917721092e-1;
const double Ry_to_eV = 13.605693009;

// Ry to He
const double Ry_to_He = 0.5;
} //end of namespace physics

namespace math{

} //end of namespace math

} //end of namespace celerium

#endif
