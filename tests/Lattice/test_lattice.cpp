#include <iomanip>      // std::setprecision
#include <lattice.h>
#include <potential.h>

using namespace celerium;

int main(int argc, char *argv[])
{

  // Load pseudopotential.
  LocalPotential potential_cr("../Potential/Cr.UPF");

  // Initialize element using chromium pseudopotential.
  Element chromium("Cr", potential_cr);

  // Set up the solver.
  wf_solver_params params;
  params.r_min = 1e-6;
  params.energy_step = 0.01;
  params.matrix_dim = 2000;
  params.grid_size = 10000;
  params.matching_index = 5000;
  params.energy_accuracy = 1e-8;

  // Add 3d-like orbitals. Note that mesh is not passed as an argument here
  // as it is provided by LocalPotential object.
  chromium.AddOrbitalClass(3,         // n = 3
                           2,         // l = 2
                           params);   // solver parameters

  // Add 4s-like orbitals.
  chromium.AddOrbitalClass(4,         // n = 4
                           0,         // l = 0
                           params);   // solver parameters

  ElementaryCell elementary_cell;

  elementary_cell.AddSite("Ce(1)", chromium, {{0, 0, 0}});
  


  std::vector<double> result;

  ArithmeticVector v({20.1, 0.2, 0.4});
  elementary_cell.EvaluateOrbitals(v, result);

  for (size_t i  = 0; i < elementary_cell.NOrbitals(); ++i) {
    std::cout << result[i] << "\n";
  }

  elementary_cell.EvaluatePotentials(v, result);


  
  std::cout << "\n\n";
  
  for (size_t i  = 0; i < elementary_cell.NSites(); ++i) {
    std::cout << result[i] << "\n";
  }
  
  return 0;
}
