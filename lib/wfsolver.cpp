#include "../headers/wfsolver.h"

void Normalize(eigenstate_struct &eigenstate) {
  double scaling_factor = 0.0;


  for (size_t i = 0; i+1 < eigenstate.wave_function.size(); ++i) {
    double r = eigenstate.wave_function[i].x;
    double dr = eigenstate.wave_function[i+1].x - eigenstate.wave_function[i].x;
    double val = eigenstate.wave_function[i].y;
    scaling_factor += r*dr*val*val;
  }
  
  for (auto &sample : eigenstate.wave_function) {
    sample.y /= std::sqrt(scaling_factor * sample.x);
  }
  
}

double CalculateY2(double y0, double y1,
                   double g0, double g1, double g2,
                   double dx2) {
  
  const double nominator =
      2.0*y1*(1.0 - 5.0/12.0*g1*dx2) - y0*(1.0 + g0*dx2/12.0);
  const double denominator =
      1.0 + g2*dx2/12.0;
  return nominator/denominator;
}



int SolveRadialSchroedingerEqn(potential_t potential,
                               size_t l,
                               double r_min,
                               double r_max,
                               double initial_energy,
                               double energy_step,
                               size_t grid_size,
                               size_t matching_index,
                               double matching_accuracy,
                               eigenstate_struct &eigenstate) {


  if (r_min < 0)
    throw std::invalid_argument(
        "SolveRadialSchroedingerEqn: argument r_min must be positive");
  if (r_max < r_min)
    throw std::invalid_argument(
        "SolveRadialSchroedingerEqn: argument \
r_max must be greater than r_min");
  if (grid_size < 10)
    throw std::invalid_argument(
        "SolveRadialSchroedingerEqn: argument grid_size \
must be grater than 9.");
  if (matching_index < 3)
    throw std::invalid_argument(
        "SolveRadialSchroedingerEqn: argument matching_index \
must be grater than 2.");
  if (matching_index + 3 > grid_size)
    throw std::invalid_argument(
        "SolveRadialSchroedingerEqn: argument matching_index must \
be smaller than grid_size - 3.");
  
  const double x_min = std::log(r_min);
  const double x_max = std::log(r_max);  
  const double dx2 = pow((x_max - x_min) / (grid_size-1.0), 2);

  eigenstate.wave_function.resize(grid_size);

  // Define grid of points.
  for (size_t i = 0; i < grid_size; ++i) {
    eigenstate.wave_function[i].x = exp(x_min+(x_max-x_min)*i/(grid_size-1.0));
  }  

  size_t n_nodes = 0;

  double energy0 = initial_energy;
  double energy1 = initial_energy + energy_step;
  double energy = 0.0;
  double val0 = 0.0;
  double val1 = 0.0;
      
  for (size_t iter = 0; iter < 200; ++iter) {

    n_nodes = 0;


    if (iter == 0) {
      energy  = initial_energy;
    }
    else if (iter == 1) {
      energy = initial_energy + energy_step;
    }
    else {
      energy = energy1 - val1/(val1 - val0)*(energy1 - energy0);
    }
  
    double r0, r1, r2, y0, y1;

    //if (iter  == 0) {
    eigenstate.wave_function[0].y = 0.1;
    eigenstate.wave_function[1].y =
        eigenstate.wave_function[0].y *
        std::pow(
            eigenstate.wave_function[1].x/eigenstate.wave_function[0].x,
            l + 0.5);
    //}
    
    for (size_t i = 2; i <= matching_index; ++i) {
      r0 = eigenstate.wave_function[i-2].x;
      r1 = eigenstate.wave_function[i-1].x;
      r2 = eigenstate.wave_function[i].x;
      y0 = eigenstate.wave_function[i-2].y;
      y1 = eigenstate.wave_function[i-1].y;
    
      
      eigenstate.wave_function[i].y =
          CalculateY2(y0,
                      y1,
                      0.262468426082*r0*r0*(-potential(r0)+energy)-pow(l+0.5,2),
                      0.262468426082*r1*r1*(-potential(r1)+energy)-pow(l+0.5,2),
                      0.262468426082*r2*r2*(-potential(r2)+energy)-pow(l+0.5,2),
                      dx2);
      if (eigenstate.wave_function[i].y * eigenstate.wave_function[i-1].y < 0) {
        n_nodes++;
      }
    }

    //if (iter  == 0) {
    eigenstate.wave_function[grid_size-1].y = 1e-6;
    eigenstate.wave_function[grid_size-2].y = 1e-6 + 1e-4*std::sqrt(dx2);
    //}
    
    for (size_t i = grid_size-3; i > matching_index; --i) {
      r0 = eigenstate.wave_function[i+2].x;
      r1 = eigenstate.wave_function[i+1].x;
      r2 = eigenstate.wave_function[i].x;
      y0 = eigenstate.wave_function[i+2].y;
      y1 = eigenstate.wave_function[i+1].y;
      
      eigenstate.wave_function[i].y =
          CalculateY2(y0,
                      y1,
                      0.262468426082*r0*r0*(-potential(r0)+energy)-pow(l+0.5,2),
                      0.262468426082*r1*r1*(-potential(r1)+energy)-pow(l+0.5,2),
                      0.262468426082*r2*r2*(-potential(r2)+energy)-pow(l+0.5,2),
                      dx2);
      
      if (eigenstate.wave_function[i].y * eigenstate.wave_function[i+1].y < 0) {
        n_nodes++;
      }
    }

    r0 = eigenstate.wave_function[matching_index+2].x;
    r1 = eigenstate.wave_function[matching_index+1].x;
    r2 = eigenstate.wave_function[matching_index].x;
    y0 = eigenstate.wave_function[matching_index+2].y;
    y1 = eigenstate.wave_function[matching_index+1].y;
    
    double y_up =
          CalculateY2(y0,
                      y1,
                      0.262468426082*r0*r0*(-potential(r0)+energy)-pow(l+0.5,2),
                      0.262468426082*r1*r1*(-potential(r1)+energy)-pow(l+0.5,2),
                      0.262468426082*r2*r2*(-potential(r2)+energy)-pow(l+0.5,2),
                      dx2);


    double scaling_factor =  eigenstate.wave_function[matching_index].y/y_up;

    for (size_t i = matching_index+1; i < grid_size; ++i) {
      eigenstate.wave_function[i].y *= scaling_factor;
    }

    Normalize(eigenstate);

    double derivative_jump =
        eigenstate.wave_function[matching_index-1].y +
        eigenstate.wave_function[matching_index+1].y -
        eigenstate.wave_function[matching_index].y *
        (2.0-(0.262468426082*r2*r2*(-potential(r2)+energy)-pow(l+0.5,2))*dx2);

    derivative_jump /= std::sqrt(dx2);
    

    
    val0 = val1;
    val1 = derivative_jump;

    energy0 = energy1;
    energy1 = energy;


    double estimated_error =
        fabs(val1*(energy - energy0)/(derivative_jump - val0));
    
    if (iter > 1 && estimated_error < matching_accuracy) {
      if  (fabs(derivative_jump) > 0.01) return -1;
      eigenstate.energy = energy;
      eigenstate.nodes = n_nodes;
      eigenstate.energy_error = fabs(energy0 - energy1);
      return 0;
    }
    
  }
  
  return -1;
}


void FindEigenvalueCandidates(
    potential_t potential,
    size_t l,
    double r_min,
    double r_max,
    size_t matrix_dim,
    std::vector<double> &eigenvalue_candidates) {

  eigenvalue_candidates.clear();

  const double dr2 = pow((r_max - r_min) / (matrix_dim-1.0), 2);
  
  gsl_matrix *m = gsl_matrix_alloc(matrix_dim, matrix_dim);
  gsl_matrix_set_zero(m);
  
  for (size_t i = 0; i < matrix_dim; ++i) {
    double r = (r_min+(r_max-r_min)*i/(matrix_dim-1.0));

    if (i + 1 < matrix_dim)
      gsl_matrix_set(m, i+1, i, -1.0/dr2/0.262468426082);
    if (i > 0)
      gsl_matrix_set(m, i-1, i, -1.0/dr2/0.262468426082);
   
    gsl_matrix_set(m, i, i, l*(l+1)/0.262468426082/r/r
                   + potential(r)
                     + 2.0/dr2/0.262468426082);   
  }  

  gsl_vector *eval = gsl_vector_alloc (matrix_dim);

  gsl_eigen_symm_workspace * w = 
    gsl_eigen_symm_alloc (matrix_dim);
  
  gsl_eigen_symm (m, eval, w);

  gsl_eigen_symm_free (w);



    for (size_t i = 0; i < matrix_dim; i++)
      {
        double eval_i 
           = gsl_vector_get (eval, i);
        if (eval_i < 0)
          eigenvalue_candidates.push_back(eval_i);
      }
  

  gsl_vector_free (eval);
  gsl_matrix_free (m);

  std::sort(eigenvalue_candidates.begin(), eigenvalue_candidates.end());
}

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
                   std::string &log) {

  std::stringstream sslog;

    if (n == 0)
    throw std::invalid_argument(
        "FindEigenstate: Invalid level index n = 0. The smallest \
allowed value is n = 1.");
    if (l+1 > n)
    throw std::invalid_argument(
        "FindEigenstate: Invalid index l. It is required that l < n.");

  
  
  std::vector<double> eigenvalue_candidates;

  sslog << "\nDetermination of candidates for bound-state eigenvalues...";
  
  FindEigenvalueCandidates(
      potential,
      l,
      r_min,
      r_max,
      matrix_dim,
      eigenvalue_candidates);

    sslog << "\nCandidates found: ";

    for (size_t i = 0; i < eigenvalue_candidates.size(); ++i) {
      sslog << eigenvalue_candidates[i];
      if (i + 1 != eigenvalue_candidates.size()) sslog << ", ";
    }

  int status;
  for(size_t candidate_number = 0;
      candidate_number < eigenvalue_candidates.size();
      ++candidate_number) {
    
    status = SolveRadialSchroedingerEqn(potential,
                                        l,
                                        r_min,
                                        r_max,
                                        eigenvalue_candidates[candidate_number],
                                        energy_step,
                                        grid_size,
                                        matching_index,
                                        eigenvalue_accuracy,
                                        eigenstate);

    if (status != 0) {
      sslog << "\nChecking eigenvalue candidate "
            << candidate_number << ":\n";
      sslog << "FAILED TO CONVERGE";
      continue;
    }

    if (eigenstate.nodes <= n - l - 1) {
      sslog << "\nChecking eigenvalue candidate "
            << candidate_number << ":\n";
      sslog << "CONVERGED (intial value: "
            << eigenvalue_candidates[candidate_number] <<", ";
      sslog << "final energy = " << eigenstate.energy
            << " += " << eigenstate.energy_error<< ", ";
      sslog << "nodes = " << eigenstate.nodes << ")";
    }
    
    if (eigenstate.nodes > n - l - 1) {
      sslog << "\n*** FAIL! Target wave function could not be found... ***";
      log = sslog.str();
      return -1;
    }

    if (eigenstate.nodes == n - l - 1) {
            sslog << "\n*** SUCCESS! Target eigenstate found... ***";
      log = sslog.str();
      return 0;
    }
  }

  sslog << "\n*** FAIL! Target wave function could not be found... ***";
  log = sslog.str();
  return -1;
}


