#include <lattice.h>
#include <potential.h>

using namespace celerium;

int main(int argc, char *argv[])
{  
  // Load pseudopotential.
  LocalPotential potential_cr("../Potential/Cr.UPF");

  // Initialize element using chromium pseudopotential.
  Element chromium("Cr",            // Element name.
                   potential_cr,    // Potential.
                   {});             // Do not add any orbital classes manually.
                                    // We will obtain orbitals numerically.

  // Now add 3d-like orbitals. Note that sampling mesh is not passed as an 
  // argument here as it is provided by LocalPotential object. Also, solver 
  // parameters are not specified (default values are used).
  chromium.AddOrbitalClass(3,         // n = 3
                           2);        // l = 2

  // Add 4s-like orbitals.
  chromium.AddOrbitalClass(4,         // n = 4
                           0);        // l = 0
  
  // Define basis vectors for Cr.
  double a {2.91};
  std::array<ArithmeticVector, 3> basis;
  basis[0] = {{a, 0.0, 0.0}};
  basis[1] = {{0.0, a, 0.0}};
  basis[2] = {{0.0, 0.0, a}};

  ElementaryCell elementary_cell(basis);
  
  elementary_cell.AddSite("Cr(1)", chromium, {{0, 0, 0}});
  elementary_cell.AddSite("Cr(2)", chromium, {{a*0.5, a*0.5, a*0.5}});

  Lattice lattice(elementary_cell);

  cuba::Cuba engine(1e9,1e6,1e-4);
  
  // Uncomment the below region to recalculate Wannier coefficients for Cr.
  /*
  std::vector<std::pair<double, double>> integration_limits(3, {-15, 15});

  lattice.CalculateWannierData(
      {{2, 2, 2}},            // Range of orbital overlaps (in the units of lattice parameters).
      {{6, 6, 6}},            // Maximum range of orbitals contributing to the Wannier functions.
      {{0, 0, 0}, {1, 0, 0}}, // Wannier located at points (0, 0, 0) and (1, 0, 0)) will be caluclated.
      integration_limits,     // Intergation limits for orbital overlaps.
      engine,                 // Integration engine.
      0.01,                   // Lower cutoff for the Wannier coefficients.
      true,                   // Suppress coefficiens below numerical accuracy treshold.
      true);                  // Verbose = true. Pass info about progress to std::cerr
                              // to allow for progress monitoring.
  
  lattice.SaveWannierDataToFile("wanniers.dat"); // Save Wannier data to file for later use.
  return 0;
  */
  
  lattice.LoadWannierDataFromFile("wanniers.dat");

  size_t npt = 100;
   for (size_t i  = 1; i < npt; ++i) {
     for (size_t j  = 1; j < npt; ++j) {
      
            double x = -5.012 + 10.0/npt*i;
            double y = -5.012 + 10.0/npt*j;

            // Update all Wanniers. This might be a time-consuming operation.
            // Once it is done, we will have a rapid access to Wannier
            // functions.
            lattice.UpdateWanniers({{x, y, 0.0}});

            // Print Wannier #5 (4s-like on Cr(1)), located at
            // position #0 (R = [0, 0, 0]).
            std::cout << x << " " << y << " " <<
                lattice.GetWannier(0, 5) << "\n";
    }
  }

  return 0;
}
