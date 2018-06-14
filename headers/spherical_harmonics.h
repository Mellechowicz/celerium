#include <complex>
#include <gsl/gsl_sf_legendre.h>

// Provides complex and real spherical harmonics.

#ifndef SPHERICAL_HARMONICS_H
#define SPHERICAL_HARMONICS_H

namespace celerium {

// The Condon–Shortley phase is NOT included
// in the spherical harmonic definition, i.e.
// Y_l^{-m}(theta, phi) = (-1)^m * Y_l^m(theta, phi)^\dagger.
std::complex<double> SphericalHarmonic(int l, int m,
                                       double theta, double phi);

// Calculates real sperical harmonic as defined here:
// https://en.wikipedia.org/wiki/Spherical_harmonics#Real_form
double RealSphericalHarmonic(int l, int m,
                             double theta, double phi);

// Calculates real sperical harmonic
// taking cartesian coordinates as input.
// Angles theta and phi are automatically extracted 
// from the coordinates in the following way:
// x = r sin(theta) cos(phi),
// y = r sin(theta) sin(phi),
// z = r cos(theta).
// The value of the spherical harmonic
// does not depend on r = sqrt(x^2 + y^2 + z^2).
double RealSphericalHarmonic(int l, int m,
                             double x,
                             double y,
                             double z);

// Include sources.
#include "../lib/spherical_harmonics.cpp"

} // end namespace celerium

#endif /* SPHERICAL_HARMONICS_H */
