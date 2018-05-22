#ifndef WFSOLVERINTERFACE_H
#define WFSOLVERINTERFACE_H

// This is a high-level interface for WFSolver library
// that solves radial Schroedinger equation for a
// user-provided radial poential.
// To use it include this header file and
// link against $(top_srcdir)/lib/WFSolver/libwfsolver.a.
// This interace provides a template WFSolver
// and the structure wf_solver_params used to
// set up the solver. For exemplary uses
// consult WFSolver in tests directory.

#include "../lib/WFSolver/include/wfsolver.h"

namespace celerium{

// Description of the parameters:
// r_min:           lower cutoff for the radius. This is required due to 
//                  logarithmic grid used by the solver.
//                  Typical values: ~ 1e-3 - 1e-12.
// energy_step:     initial energy step in the eigenvalue search algorithm
//                  typical values: 0.001-0.1
// matrix_dim:      dimension of a matrix eigenvalue problem used to roughly
//                  estimate the eigenvalues.
//                  Typical values: 1000-2000.
// grid_size:       size of a mesh used for finding accurate solutions of
//                  Schroedinger equation.
//                  Typical values: 10^3 - 10^5.
// matching_index:  index for which two iterative procedures starting from two
//                  sides of the system are matched.
//                  Typical values: around grid_size/2.
// energy_accuracy: stopping condition for eigenenergy determiantion.
//                  The procedure stops if the enegy difference between
//                  two consecutive steps is smaller than energy_accuracy.
//                  In the case of slow convergence, energy_accuracy
//                  might not reflect actual accuracy of computed
//                  eigenvalue.
//                  Typical values: 1e-6 - 1e-10
struct wf_solver_params {
  double r_min;
  double energy_step;
  size_t matrix_dim;
  size_t grid_size;
  size_t matching_index;
  double energy_accuracy;
};


template <class Potential>
class WFSolver {
  
 public:

  // Initializes the solver with a potential
  // and parameters.
  WFSolver() {
    this->logs = "";
  }

  // Return logs from the last call of GetEigenstate method.
  // These should be consulted if GetEigenstate returns
  // a non-zero value.
  const std::string &GetLog() const {return this->logs;}

  // Calculates wave function and energy for the
  // n-th eigenstate with angular momentum l.
  // In case of success returns zero and
  // a non-zero value otherwise.
  // The eigenstate numbering starts from ONE, i.e.
  // n = 1, 2, 3 ..., whereas l = 0, 1, 2, ...
  // for the s, p, d, ... states.
  int GetEigenstate(size_t n, size_t l,
                    Potential &potential,
                    const wf_solver_params &params,
                    std::vector<double> &wave_function,
                    double &energy) {

    wave_function.clear();

    std::vector<double> mesh;
    std::vector<double> values;
    potential.get_mesh(mesh);
    potential.get_local(values);

    std::vector<sample_struct> samples(mesh.size());
    for (size_t i  = 0; i < mesh.size(); ++i) {
      samples[i].x = mesh[i];
      samples[i].y = values[i];
    }
       
    double r_max = mesh.front();
    for (auto r : mesh) {
      if (r > r_max) r_max = r;
    }

    Interpolator interp(samples, gsl_interp_cspline);

    eigenstate_struct eigenstate;
    int status = FindEigenstate(interp,
                                n,
                                l,
                                params.r_min,
                                r_max,
                                params.matrix_dim,
                                params.energy_step,
                                params.grid_size,
                                params.matching_index,
                                params.energy_accuracy,
                                eigenstate,
                                this->logs);
    
      if (status == 0) {
        energy = eigenstate.energy;
        Interpolator wave_function_interp(eigenstate.wave_function,
                                          gsl_interp_cspline);
        for (double r : mesh) {
          wave_function.push_back(wave_function_interp(r));
        }      
        return 0;
      }
    
    return -1;
  }

 private:
  std::string logs;
};

}

#endif /* WFSOLVERINTERFACE_H */

