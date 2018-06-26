#ifndef ELEMENT_H
#define ELEMENT_H

#include <vector>
#include <string>
#include "wfsolver_interface.h"
#include "potential.h"
#include "orbital_class.h"

namespace celerium {

template <class RadialPotential = Interpolator,
          class OrbitalClass = OrbitalClass<Interpolator>> 
class Element {
 public:

  Element(const char* name,
          const RadialPotential &radial_potential,
          const std::vector<OrbitalClass> &orbital_classes) :
      radial_potential(radial_potential),
      name(name),
      orbital_classes(orbital_classes) {}

  Element(const char* name,
          RadialPotential &&radial_potential,
          const std::vector<OrbitalClass> &orbital_classes) :
      radial_potential(std::move(radial_potential)),
      name(name),
      orbital_classes(orbital_classes) {}

  Element(const char* name,
          const RadialPotential &radial_potential,
          const std::initializer_list<OrbitalClass> &orbital_classes) :
      radial_potential(radial_potential),
      name(name),
      orbital_classes(orbital_classes) {}

  Element(const char* name,
          RadialPotential &&radial_potential,
          const std::initializer_list<OrbitalClass> &orbital_classes) :
      radial_potential(std::move(radial_potential)),
      name(name),
      orbital_classes(orbital_classes) {}

  
  Element(const char *element_name,
          LocalPotential &potential) {
    this->name = element_name;
    this->SetRadialPotential(potential);
  }
  
  void SetRadialPotential(LocalPotential &potential) {
    std::vector<sample_struct> samples;
    std::vector<double> x, y;
    potential.get_mesh(x);
    potential.get_local(y);   
    for (size_t i  = 0; i < x.size(); ++i) samples.push_back({x[i], y[i]});
    this->radial_potential = Interpolator(samples);
  }
    
  int AddOrbitalClass(int n, int l,
                      const std::vector<int> &active_m_values,
                      const wf_solver_params &params) {
    
    double r_max = this->radial_potential.GetSamples().front().x;
    for (auto s : this->radial_potential.GetSamples()) {
      if (s.x > r_max) r_max = s.x;
    }                  
    
    std::string logs;
    eigenstate_struct eigenstate;
    int status = FindEigenstate(this->radial_potential, n, l,
                                params.r_min, r_max,
                                params.matrix_dim,
                                params.energy_step,
                                params.grid_size,
                                params.matching_index,
                                params.energy_accuracy,
                                eigenstate,
                                logs);
    
    if (status == 0) {

      Interpolator interp_result(eigenstate.wave_function);
      std::vector<sample_struct> new_samples;
      
      for (const auto &s : this->radial_potential.GetSamples()) {
        new_samples.push_back({s.x, interp_result(s.x)});
      }

      Interpolator wave_function_original_mesh_interp(new_samples);
      OrbitalClass orbital_class(wave_function_original_mesh_interp,
                                 eigenstate.energy,
                                 l, active_m_values);

      this->orbital_classes.push_back(orbital_class);
      return 0;
    }
    else {

      std::cout << "\nSOLVER FAILED FOR ELEMENT"
                << this->name << ", l = " << l << ".\n";
      std::cout << "\n\n" << logs << "\n\n";  
      throw std::runtime_error("celerium::Element::AddOrbital: Schroedinger \
equation solver did not converge (adjust the solver settings).");
    }
    return -1;
  }

  void AddOrbitalClass(int n, int l, const wf_solver_params &params) {
    std::vector<int> active_m_values;
    for (int m = -l; m <= l; ++m) active_m_values.push_back(m);
    this->AddOrbitalClass(n, l, active_m_values, params);
  }
  
  // Setters.

  void SetRadialPotential(RadialPotential &radial_potential) {
    this->radial_potential = radial_potential;
  }

  void SetName(const char *element_name) {
    this->name = element_name;
  }

  void SetOrbitalClasses(std::vector<OrbitalClass> &orbital_classes) {
    this->orbital_classes = orbital_classes;
  }

  // Getters.
  
  const std::string &GetName() const {return this->name;}

  RadialPotential &GetRadialPotential(){return this->radial_potential;}

  std::vector<OrbitalClass> &GetOrbitalClasses() {
    return this->orbital_classes;
  }

 private:

  std::string name;
  RadialPotential radial_potential;
  std::vector<OrbitalClass> orbital_classes;
  
};





}  // end namespace celerium







/*

template <class RadialWaveFunction>
struct extended_radial_wf {
  RadialWaveFunction radial_wf;
  int l;
  double energy;
};



template <class RadialPotential, class RadialWaveFunction>
class Element {

 public:

  Element() {};
  
  Element(const std::string &name, const RadialPotential &radial_potential) {
    this->name = name;
    this->radial_potential = radial_potential;
  }
  
  void SetName(const std::string &name) {this->name = name;}

  void SetRadialPotential(const RadialPotential &radial_potential) {
    this->radial_potential = radial_potential;
  }

  void AddRadialWaveFunction(RadialWaveFunction radial_wave_function,
                             int l,
                             double energy) {
    
    if (l < 0) throw std::invalid_argument("celerium::Element::AddRadialWaveFunction:\
 l must by greater of equal to zero.");

    auto located_wf = std::find_if(this->extended_radial_wfs.begin(),
                                   this->extended_radial_wfs.end(),
                                   [&](extended_radial_wf<RadialWaveFunction>
                                       tested_extended_wf) {
                                     return tested_extended_wf.l == l;
                                   });

    if (located_wf != this->extended_radial_wfs.end())
      throw std::invalid_argument("celerium::Element::AddWaveFunction:\
 attempted to add two wave functions corresponding to the same l.");

    this->extended_radial_wfs.push_back({radial_wave_function, l, energy});
  }

  double GetWaveFunction(ArithmeticVector &coords,
                         size_t radial_wf_index,
                         int m) {

    auto &extended_radial_wf = this->extended_radial_wfs[radial_wf_index];
    
    if (abs(m) > extended_radial_wf.l)
      throw std::invalid_argument("celerium::Element::GetWaveFunction\
 |m| must not exceed l.");    

    double result = extended_radial_wf.radial_wf(coords.length());
    result *= RealSphericalHarmonic(extended_radial_wf.l, m,
                                    coords[0],
                                    coords[1],
                                    coords[2]);
    return result;
  }

  const std::string &GetName() const {return this->name;}

  double GetRadialPotential(double r) {
    return this->radial_potential(r);
  }

  const std::vector<extended_radial_wf<RadialWaveFunction>> &
  GetRadialWaveFunctions() const {
    return this->extended_radial_wfs;
  }
  
 private:
  std::string name;
  RadialPotential radial_potential;
  std::vector<extended_radial_wf<RadialWaveFunction>> extended_radial_wfs;
};


using ElementI = Element<Interpolator, Interpolator>;


using ElementF = Element<std::function<double(double)>,
                         std::function<double(double)>>;




} // end namespace celerium

*/

#endif /* ELEMENT_H */

