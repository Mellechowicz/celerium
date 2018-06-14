#include "../../headers/potential.h"
#include "../../headers/wfsolver_interface.h"
#include <iomanip>

using namespace celerium;

int main() {

  LocalPotential potential;
  potential.input("../Potential/Cr.UPF");

  WFSolver<LocalPotential> solver;

  double energy;
  std::vector<double> wave_function;

  wf_solver_params params;  
  params.r_min = 1e-7;
  params.energy_step = 0.01;
  params.matrix_dim = 2000;
  params.grid_size = 5000;
  params.matching_index = 2500;
  params.energy_accuracy = 1e-8;
  
  int status = solver.GetEigenstate(3, 1,
                                    potential,
                                    params,
                                    wave_function,
                                    energy);


  std::cout << "\nDetermination of the 3p-like state for Cr:";
  std::cout << "\n\n" << solver.GetLog() << "\n\n";
  
  if (status == 0) {

    std::cout << std::fixed << std::setprecision(12);
    std::cout << "\nCalculated energy: " << energy << " eV\n\n";
    
    std::vector<double> potential_values;
    std::vector<double> mesh;

    potential.get_local(potential_values);
    potential.get_mesh(mesh);

    std::cout << "r               V(r)            Psi(r) \n";
    
    for (size_t i = 0; i < mesh.size(); ++i) {
      std::cout << mesh[i] << " " << potential_values[i]
                << " " << wave_function[i] << "\n";
    }

  }
  else {
    std::cout << "FAILED TO CONVERGE!";
    
    return -1;
  }
  
  return 0;
}
