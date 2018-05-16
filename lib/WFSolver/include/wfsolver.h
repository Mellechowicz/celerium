#ifndef WFSOLVER_H
#define WFSOLVER_H

#include <functional>
#include <vector>
#include <iostream>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_eigen.h>
#include <algorithm>
#include "interpolator.h"
#include <sstream>


typedef std::function<double(double)> potential_t;

struct eigenstate_struct {
  double energy;
  double energy_error;
  size_t nodes;
  std::vector<sample_struct> wave_function;
};

int SolveRadialSchroedingerEqn(potential_t potential,
                                size_t l,
                                double r_min,
                                double r_max,
                                double initial_energy,
                                double energy_step,
                                size_t grid_size,
                                size_t matching_index,
                                double matching_accuracy,
                                eigenstate_struct &eigenstate);

void FindEigenvalueCandidates(
    potential_t potential,
    size_t l,
    double r_min,
    double r_max,
    size_t matrix_dim,
    std::vector<double> &eigenvalue_candidates);



int FindEigenstate(potential_t potential,
                   size_t n,
                   size_t l,
                   double r_min,
                   double r_max,
                   size_t matrix_dim,
                   double energy_step,
                   size_t grid_size,
                   size_t matching_index,
                   double eigenvalue_accuracy,
                   eigenstate_struct &eigenstate,
                   std::string &log);





#endif /* WFSOLVER_H */

