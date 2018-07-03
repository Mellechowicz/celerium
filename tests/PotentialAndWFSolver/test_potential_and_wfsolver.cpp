#include "../../headers/potential.h"
#include "../../headers/wfsolver.h"
#include <iomanip>

using namespace celerium;

int main() {

  LocalPotential potential;
  potential.input("../Potential/Cr.UPF");

  WFSolver solver;

  double energy;

  eigenstate_struct eigenstate;
  std::vector<double> mesh, values;

  potential.get_mesh(mesh);
  potential.get_local(values);
  
  int status = solver.GetEigenstate(3, 1, mesh, values, mesh, eigenstate);

  std::cout << "\nDetermination of the 3p-like state for Cr:";
  std::cout << "\n\n" << solver.GetLogs() << "\n\n";
  
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
                << " " << eigenstate.wave_function[i].y << "\n";
    }

  }
  else {
    std::cout << "FAILED TO CONVERGE!";
    
    return -1;
  }
  
  return 0;
}

