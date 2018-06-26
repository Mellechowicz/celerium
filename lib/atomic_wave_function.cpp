#include <atomic_wave_function.h>

template<class RadialWaveFunction>
AtomicWaveFunction<RadialWaveFunction>::AtomicWaveFunction()
{
  this->l = 0;
  this->m = 0;
}

template<class RadialWaveFunction>
AtomicWaveFunction<RadialWaveFunction>::AtomicWaveFunction(
    RadialWaveFunction radial_function,
    int l,
    int m) {

  if (abs(m) > l)
    throw std::invalid_argument("AtomicWaveFunction: m cannot exceed l.");
  
  this->radial_function = radial_function;
  this->l = l;
  this->m = m;
}

template<class RadialWaveFunction>
void AtomicWaveFunction<RadialWaveFunction>::SetRadialFunction(
    RadialWaveFunction radial_function) {
  this->radial_function = std::move(radial_function);
}

template<class RadialWaveFunction>
void AtomicWaveFunction<RadialWaveFunction>::SetQuantumNumbers(
    size_t l, size_t m) {

  if (abs(m) > l)
    throw std::invalid_argument("AtomicWaveFunction.SetQuantumNumbers: m cannot exceed l.");
  
  this->l = l;
  this->m = m;
}

template<class RadialWaveFunction>
double AtomicWaveFunction<RadialWaveFunction>::operator()(
    const ArithmeticVector &cartesian_coords) {

  const double x = cartesian_coords[0];
  const double y = cartesian_coords[1];
  const double z = cartesian_coords[2];
  const double r = sqrt(cartesian_coords[0]*cartesian_coords[0] +
                  cartesian_coords[1]*cartesian_coords[1] +
                  cartesian_coords[2]*cartesian_coords[2]);
  
  double result = this->radial_function(r);
  result *= RealSphericalHarmonic(this->l, this->m, x, y, z);
 
  return  result;
}

template<class RadialWaveFunction>
double AtomicWaveFunction<RadialWaveFunction>::Laplacian(
    const ArithmeticVector &cartesian_coords) {


  const double r = sqrt(cartesian_coords[0]*cartesian_coords[0] +
                  cartesian_coords[1]*cartesian_coords[1] +
                  cartesian_coords[2]*cartesian_coords[2]);

  if (r == 0.0) return 0.0;

  double result = this->radial_function.D2(r) +
                  this->radial_function.D1(r)*2.0/r -
                  1.0/r/r * this->l*(this->l + 1) * this->radial_function(r);

  result *= RealSphericalHarmonic(this->l, this->m,
                                  cartesian_coords[0],
                                  cartesian_coords[1],
                                  cartesian_coords[2]);
  
  return result;  
}
