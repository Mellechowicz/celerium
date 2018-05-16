#ifndef WFSOLVERINTERFACE_H
#define WFSOLVERINTERFACE_H

#include "../lib/WFSolver/include/wfsolver.h"

namespace celerium{

// Description of the parameters:
// r_min:           cutoff for the radius. Potential in the range (0, r_min) is
//                  assumed to be infinite. Note that for s-type solutions r_min
//                  should be small (even of the order of 1e-10).
//                  For other wave function r_min~0.001-0.1
//                  might be sufficient.
// energy_step:     initial energy step in the eigenvalue search algorithm
//                  typical values: 0.001-0.1
// matrix_dim:      dimension of a matrix eigenvalue problem used to roghly
//                  estimate the eigenvalues.
//                  Typical values: 1000-2000.
// grid_size:       size of a mesh used for finding accurate solutions of
//                  Schroedinger equation.
//                  Typical values: 10^3 - 10^5.
// matching_index:  index for which two iterative procedures starting from two
//                  sides of the syste are matched.
//                  Typical values: around grid_size/2.
// energy_accuracy: target accuracy of energy determiantion.
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
  
  WFSolver(const Potential &potential,
                     const wf_solver_params &params) {
    this->potential = potential;
    this->params = params;
    this->log = "";
  }

  const std::string &GetLog() const {return this->log;}
  
  int GetEigenstate(size_t n, size_t l,
                      std::vector<double> &wave_function,
                    double &energy){

    wave_function.clear();

    std::vector<double> mesh;
    std::vector<double> values;
    this->potential.get_mesh(mesh);
    this->potential.get_local(values);

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
                                this->params.r_min,
                                r_max,
                                this->params.matrix_dim,
                                this->params.energy_step,
                                this->params.grid_size,
                                this->params.matching_index,
                                this->params.energy_accuracy,
                                eigenstate,
                                this->log);
    
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
  Potential potential;
  wf_solver_params params;
  std::string log;
};

}

#endif /* WFSOLVERINTERFACE_H */

