#ifndef ORBITAL_H
#define ORBITAL_H

#include <stdexcept>
#include "spherical_harmonics.h"
#include "arithmeticvector.h"

namespace celerium {

template <class RadialWaveFuction>
class Orbital {

 public:

  Orbital() = default;

  Orbital(RadialWaveFuction &radial_wf, int l, double energy) {
    
    if (l < 0) throw std::invalid_argument("celerium::Orbital: \
 Attempted to construct Orbital object with negative l.");
    
    this->radial_wf = radial_wf;
    this-> energy = energy;
    this->l = l;
  }

  void SetRadialWaveFunction(RadialWaveFuction &&radial_wf) {
    this->radial_wf = radial_wf;
  }

  void SetEnergy(double energy) {
    this->energy = energy;
  }

  void SetL(int l) {

    if (l < 0) throw std::invalid_argument("celerium::Orbital:SetL \
 Attempted to set negative l.");
    
    this->l = l;
  }

  double Eval(ArithmeticVector coords, int m) {

    if (abs(m) > this->l)
      throw std::invalid_argument("celerium::Orbital::Eval:\
 |m| must not exceed l.");    

    double result = this->radial_wf(coords.length());
    result *= RealSphericalHarmonic(this->l, m,
                                    coords[0],
                                    coords[1],
                                    coords[2]);
    return result;
  }
  
 private:
  RadialWaveFuction radial_wf;
  int l;
  double energy;
};


#endif /* ORBITAL_H */
