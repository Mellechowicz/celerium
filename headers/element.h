#ifndef ELEMENT_H
#define ELEMENT_H

#include <vector>
#include <string>
#include "wfsolver.h"
#include "potential.h"
#include "orbital_class.h"
#include "interpolator.h"

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
          LocalPotential &potential,
          const std::initializer_list<OrbitalClass> &orbital_classes) {
    this->name = element_name;
    this->SetRadialPotential(potential);
    this->orbital_classes = orbital_classes;
  }
  
  void SetRadialPotential(LocalPotential &potential) {
    std::vector<sample_struct> samples;
    std::vector<double> x, y;
    potential.get_mesh(x);
    potential.get_local(y);   
    for (size_t i  = 0; i < x.size(); ++i) samples.push_back({x[i], y[i]});
    this->radial_potential = Interpolator(samples);
  }
  
  void AddOrbitalClass(int n, int l,
                       const std::vector<int> &active_m_values,
                       const wf_solver_params_struct &params,
                       const std::vector<double> &result_mesh) {

    if (result_mesh.size() == 0)
      throw std::invalid_argument("celerium::Element::AddOrbitalClass: \
mesh cannot be empty.");

    double r_max = result_mesh.front();
    for (auto r : result_mesh) {

      if (r < 0)
        throw std::invalid_argument("celerium::Element::AddOrbitalClass: \
mesh may contain only non-negative numbers.");
      
        if (r > r_max) r_max = r;
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
      std::vector<sample_struct> mesh_samples;
      
      for (auto r : result_mesh) 
        mesh_samples.push_back({r, interp_result(r)});

      Interpolator wave_function_original_mesh_interp(mesh_samples);
      OrbitalClass orbital_class(wave_function_original_mesh_interp,
                                 eigenstate.energy,
                                 n, l, active_m_values);

      this->orbital_classes.push_back(orbital_class);
    }
    else {

      std::cout << "\nSOLVER FAILED FOR ELEMENT"
                << this->name << ", l = " << l << ".\n";
      std::cout << "\n\n" << logs << "\n\n";  
      throw std::runtime_error("celerium::Element::AddOrbital: Schroedinger \
equation solver did not converge (adjust the solver settings).");
    }
 
  }

  void AddOrbitalClass(int n, int l,
                       const std::vector<int> &active_m_values,
                       const std::vector<double> &result_mesh) {

    wf_solver_params_struct params;
    AddOrbitalClass(n, l, active_m_values, params, result_mesh);
  }
  
  void AddOrbitalClass(int n, int l,
                       const wf_solver_params_struct &params,
                       const std::vector<double> &result_mesh) {
    std::vector<int> active_m_values;
    for (int m = -l; m <= l; ++m) active_m_values.push_back(m);
    this->AddOrbitalClass(n, l, active_m_values, params, result_mesh);
  }

  void AddOrbitalClass(int n, int l,
                       const std::vector<double> &result_mesh) {
    std::vector<int> active_m_values;
    for (int m = -l; m <= l; ++m) active_m_values.push_back(m);
    this->AddOrbitalClass(n, l, active_m_values, result_mesh);
  }

  void AddOrbitalClass(int n, int l,
                      const std::vector<int> &active_m_values,
                      const wf_solver_params_struct &params) {

    std::vector<sample_struct> samples = this->radial_potential.GetSamples();
    std::vector<double> mesh;
    for (const auto &s : samples) mesh.push_back(s.x);
    this->AddOrbitalClass(n, l, active_m_values, params, mesh);
  }

  void AddOrbitalClass(int n, int l,
                      const std::vector<int> &active_m_values) {

    wf_solver_params_struct params;
    this->AddOrbitalClass(n, l, active_m_values, params);
  }


  void AddOrbitalClass(int n, int l, const wf_solver_params_struct &params) {
    std::vector<int> active_m_values;
    for (int m = -l; m <= l; ++m) active_m_values.push_back(m);
    this->AddOrbitalClass(n, l, active_m_values, params);
  }

  void AddOrbitalClass(int n, int l) {
    wf_solver_params_struct params;
    this->AddOrbitalClass(n, l, params);
  }
  
  // Setters.

  void SetRadialPotential(const RadialPotential &radial_potential) {
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

  const RadialPotential &GetRadialPotential() const {
    return this->radial_potential;
  }

  const std::vector<OrbitalClass> &GetOrbitalClasses() const {
    return this->orbital_classes;
  }

 private:

  RadialPotential radial_potential;
  std::string name;
  std::vector<OrbitalClass> orbital_classes;
  
};


}  // end namespace celerium

#endif /* ELEMENT_H */
