#ifndef ORBITAL_CLASS_H
#define ORBITAL_CLASS_H

#include <stdexcept>
#include "spherical_harmonics.h"
#include "arithmeticvector.h"
#include <vector>
#include <algorithm>

namespace celerium {


template <class RadialWF>
class OrbitalClass {

 public:

  OrbitalClass() = default;

  // Initializes an orbital class with guven radial wave function,
  // eigenenergy and angular orbtial momentum. All values
  // m = -l, .. l are considered active bu default.
  OrbitalClass(const RadialWF &radial_wf, double energy, int l) {
    
    if (l < 0) throw std::invalid_argument("celerium::OrbitalClass: \
 Attempted to construct Orbital object with negative l.");
    
    this->radial_wf = radial_wf;
    this-> energy = energy;
    this->l = l;

    for (int m = -l; m <= l; ++m) this->active_m_values.push_back(m);
  }

  // Initializes an orbital class with given radial wave function,
  // eigenenergy and angular orbtial momentum. Only
  // provided m values are considered active;
  OrbitalClass(const RadialWF &radial_wf, double energy,
               int l, std::vector<int> active_m_values) {
    
    if (l < 0) throw std::invalid_argument("celerium::OrbitalClass: \
 Attempted to construct Orbital object with negative l.");
    
    this->radial_wf = radial_wf;
    this-> energy = energy;
    this->l = l;
    this->SetActiveMValues(active_m_values);
  }

  // Setters.
  
  void SetRadialWaveFunction(const RadialWF &radial_wf) {
    this->radial_wf = radial_wf;
  }

  void SetEnergy(double energy) {this->energy = energy;}

  void SetL(int l) {

    if (l < 0) throw std::invalid_argument("celerium::OrbitalClass::SetL \
 Attempted to set negative l.");
    
    this->l = l;
  }
  
  // Sets active m values. Autmatically checks if |m| < l and
  // removes duplicates.
  void SetActiveMValues(const std::vector<int> &active_m_values) {
    for (auto &m : active_m_values) {

      if (abs(m) > l) throw std::invalid_argument("celerium::OrbitalClass::\
SetActiveMValues Active m values must satisfy |m| < l.");

      if ( std::find(this->active_m_values.begin(),
                     this->active_m_values.end(), m) !=
           this->active_m_values.end() ) continue;

      this->active_m_values.push_back(m);
    }
  }
  
  // Getters.
  
  int GetL() const {return this->l;}

  double GetEnergy() const {return this->energy;}

  RadialWF &GetRadialWF() {return this->radial_wf;}

  const std::vector<int> &GetActiveMValues() const {
    return this->active_m_values;
  }

  // Evaluates the the m-th orbital of the orbital class at the point coords.
  // For instance, if l=1 then m=-1, 0, and 1 correspond to
  // p_y, p_z, and p_x, respecitvely. Note that REAL spherical harmonics are
  // used in computation of the orbtials, see
  // https://en.wikipedia.org/wiki/Spherical_harmonics#Real_form
  double Eval(ArithmeticVector coords, int m) {

    if (abs(m) > this->l)
      throw std::invalid_argument("celerium::OrbitalClass::Eval:\
 |m| must not exceed l.");    

    
    double result = this->radial_wf(coords.length());
    result *=
        RealSphericalHarmonic(this->l, m, coords[0], coords[1], coords[2]);
    return result;
  }
  
 private:
  
  RadialWF radial_wf;
  int l;
  double energy;
  std::vector<int> active_m_values;
};

} // end namespace celerium

#endif /* ORBITAL_CLASS_H */