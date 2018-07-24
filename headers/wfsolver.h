#include <functional>
#include <vector>
#include <iostream>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_eigen.h>
#include <algorithm>
#include "interpolator.h"
#include <sstream>

// This is a low-level solver of the spherically-symmertic Schroedinger
// equation. The radial potential is provided through a callable object.
// The results are stored in the eigenstate_struct, containing
// radial wave function (represented as a vector of sample_struct,
// i.e., [x, \Psi(x)] pairs), eigenergy, and additional information:
// number of nodes in the radial wave function, and crudely approximated
// uncertainty of eigenvalue.

#ifndef WFSOLVER_H
#define WFSOLVER_H

namespace celerium {

typedef std::function<double(double)> potential_t;

struct wf_solver_params_struct {
  double r_min = 1e-6;
  double energy_step = 0.01;
  size_t matrix_dim = 3000;
  size_t grid_size = 10000;
  size_t matching_index = 5000;
  double energy_accuracy = 1e-8;
};

struct eigenstate_struct {
  double energy;
  double energy_error;
  int nodes;
  std::vector<celerium::gsl::sample_struct> wave_function;
};

template <class Potential>
int SolveRadialSchroedingerEqn(Potential potential,
                                int l,
                                double r_min,
                                double r_max,
                                double initial_energy,
                                double energy_step,
                                size_t grid_size,
                                size_t matching_index,
                                double matching_accuracy,
                                eigenstate_struct &eigenstate);

template <class Potential>
void FindEigenvalueCandidates(
    Potential potential,
    size_t l,
    double r_min,
    double r_max,
    size_t matrix_dim,
    std::vector<double> &eigenvalue_candidates);


template <class Potential>
int FindEigenstate(Potential potential,
                   int n,
                   int l,
                   double r_min,
                   double r_max,
                   size_t matrix_dim,
                   double energy_step,
                   size_t grid_size,
                   size_t matching_index,
                   double eigenvalue_accuracy,
                   eigenstate_struct &eigenstate,
                   std::string &log);


class WFSolver {
 public:

  int GetEigenstate(int n, int l,
                    const std::vector<celerium::gsl::sample_struct> &potential,
                     const std::vector<double> &result_mesh,
                     eigenstate_struct &eigenstate) {

    if (potential.size() == 0)
      throw std::invalid_argument("celerium::WFSolver::GetEigenstate: potential \
must not be empty.");
    
    double r_max = potential.front().x;
    for (auto s : potential)
      if (s.x > r_max) r_max = s.x;

    celerium::gsl::Interpolator interp(potential);
    
    int status = FindEigenstate(interp,
                                n,
                                l,
                                this->params.r_min,
                                r_max,
                                this->params.matrix_dim,
                                this->params.energy_step,
                                this->params.grid_size,
                                this->params.matching_index,
                                this->params.energy_accuracy,
                                eigenstate,
                                this->logs);

    gsl::Interpolator interp_wf(eigenstate.wave_function);

    eigenstate.wave_function.clear();
    
    for (auto r : result_mesh)
      eigenstate.wave_function.push_back({r, interp_wf(r)});
    
    return status;
  }
  

  int GetEigenstate(int n, int l,
                     const std::vector<double> &potential_mesh,
                     const std::vector<double> &potential_values,
                     const std::vector<double> &result_mesh,
                    eigenstate_struct &eigenstate) {

    std::vector<gsl::sample_struct> potential;

    for (size_t i = 0; i < potential_mesh.size(); ++i)
      potential.push_back({potential_mesh[i], potential_values[i]});
    
    return GetEigenstate(n, l, potential, result_mesh, eigenstate);
  }

  const std::string &GetLogs() {return this->logs;}

  wf_solver_params_struct params;

 private:
    std::string logs;
};

  


// Include sources
#include "../lib/wfsolver.cpp"

} // end namespace celerium

#endif /* WFSOLVER_H */


