#include "arithmeticvector.h"
#include "spherical_harmonics.h"
#include <stdexcept>

// AtomicWaveFunction is a class that represents the atomic orbital.
// It is initalized with the radial wave function, defined as a
// callable object, returning double. If the object provides also the
// methods D1 and D2 that return first- and second- derivatives of the
// radial function, laplacian of the full wave function can be comupted
// by using the Laplacian method. The ramaining necessary parameters are:
// 1) the orbital angular momentum l = 0, 1, 2, ..., and
// 2) index m = -l, .., l.
// The real version o spherical harmonics is used, so m does not correspond to
// the z-axis projection of angular momentum, but maps onto the standard
// orbital representations, e.g.,
// l = 1, m = -1 -> p_y,
// l = 1, m = 0 -> p_z,
// l = 1, m = 1 -> p_x.

#ifndef ATOMIC_WAVE_FUNCTION_H
#define ATOMIC_WAVE_FUNCTION_H

namespace celerium {

// Representation of the atomic orbital.
// RadialWaveFuction must provide a method
//    double operator()(double r);
// If the laplacian of the wave function
// will be calculated, it also must provide
// two other methods:
//    double D1(double r);
//    double D2(double r);
// that represent first- and second-derivative
// of the radial function, respectively.
// Note that, e.g., Interpolator object, obtained
// from numerical solution of Schroedinger
// equation provides all above methods.
template<class RadialWaveFunction>
class AtomicWaveFunction {

 public:

  AtomicWaveFunction();

  // l is the orbital angular momentum (l = 1, 2, 3, ...)
  // m = -l, .., l.
  AtomicWaveFunction(RadialWaveFunction radial_function,
                     int l,
                     int m);

  void SetRadialFunction(RadialWaveFunction radial_function);

  void SetQuantumNumbers(size_t l, size_t m);

  // Returns the value of the wave function at cartesian_coords.
  double operator()(const ArithmeticVector &cartesian_coords);  

  // Returns the laplacian of the wave function at cartesian_coords.
  double Laplacian(const ArithmeticVector &cartesian_coords);

 private:
  RadialWaveFunction radial_function;
  int l;
  int m;
};

// Include sources.
#include "../lib/atomic_wave_function.cpp"

} // end namespace celerium

#endif /* ATOMIC_WAVE_FUNCTION_H */
