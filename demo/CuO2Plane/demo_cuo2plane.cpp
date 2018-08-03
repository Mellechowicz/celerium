#include <lattice.h>
#include <potential.h>

using namespace celerium;

int main(int argc, char *argv[])
{  
  // Load pseudopotential.
  LocalPotential potential_o("O.UPF");
  LocalPotential potential_cu("Cu.UPF");

  Element oxygen("O",            // Element name.
                 potential_o,    // Potential.
                 {});            // Do not add any orbital classes manually.

  Element copper("Cu",            // Element name.
                 potential_cu,    // Potential.
                 {});             // Do not add any orbital classes manually.
  
  oxygen.AddOrbitalClass(2, 0);    // 2s      
  oxygen.AddOrbitalClass(2, 1);    // 2p      

  copper.AddOrbitalClass(4, 0);    // 4s
  copper.AddOrbitalClass(4, 1);    // 4p
  copper.AddOrbitalClass(3, 2);    // 3d
  
  // Define basis vectors for CuO2 plane.
  double a {3.81380};
  double c {13.22490};
  std::array<ArithmeticVector, 3> basis;
  basis[0] = {{a, 0.0, 0.0}};
  basis[1] = {{0.0, a, 0.0}};
  basis[2] = {{0.0, 0.0, 1000.0*a}};

  ElementaryCell elementary_cell(basis);
  
  elementary_cell.AddSite("Cu", copper, {{0, 0, 0}});
  elementary_cell.AddSite("O(1)", oxygen, {{a*0.5, 0.0, 0.0}});
  elementary_cell.AddSite("O(2)", oxygen, {{0.0, a*0.5, 0.0}});
  elementary_cell.AddSite("O(3)", oxygen, {{0.0, 0.0, c*0.18623}});
  elementary_cell.AddSite("O(4)", oxygen, {{0.0, 0.0, -c*0.18623}});
  
  

  std::vector<std::string> orbital_descriptions;
  elementary_cell.GetOrbitalDescriptions(orbital_descriptions);

  for (const auto &orbital_description : orbital_descriptions)
    std::cout << orbital_description << "\n";
      

  elementary_cell.SetCrystalPotentialCutoff(40.0);
  Lattice lattice(elementary_cell);
  

  
  
  // Uncomment the below region to recalculate Wannier coefficients for Cr.

  
  /*
    cuba::Cuba engine(1e9,1e6,1e-4);
    std::vector<std::pair<double, double>> integration_limits(3, {-15, 15});

    lattice.CalculateWannierData(
      {{2, 2, 0}},            // Range of orbital overlaps (in the units of lattice parameters).
      {{6, 6, 0}},            // Maximum range of orbitals contributing to the Wannier functions.
      {{0, 0, 0}, {1, 0, 0}}, // Wannier located at points (0, 0, 0) and (1, 0, 0)) will be caluclated.
      integration_limits,     // Intergation limits for orbital overlaps.
      engine,                 // Integration engine.
      0.0001,                 // Lower cutoff for the Wannier coefficients.
      true,                   // Suppress coefficiens below numerical accuracy treshold.
      true);                  // Verbose = true. Pass info about progress to std::cerr
                              // to allow for progress monitoring.
  
  lattice.SaveWannierDataToFile("wanniers.dat"); // Save Wannier data to file for later use.
  return 0;
  */


  
  lattice.LoadWannierDataFromFile("wanniers.dat");

  cuba::Cuba engine(1e7,1e6,1e-3);
  
  std::vector<std::pair<double,double>> b3(3,std::make_pair(-10,10));
  std::vector<double> resN (25), errN (25), pN (25);  
  int steps = 0;

  std::function<int(const double *, double *)> integrand;
  integrand = [&](const double *xx, double *ff) {

    ArithmeticVector coords({xx[0], xx[1], xx[2]});
    lattice.UpdateWanniers(coords);
    lattice.UpdateLaplacians(coords);

    double w, dw, v;

    v = lattice.EvaluateCrystalPotential(xx);
    
    for (size_t i = 0; i < 24; ++i) {
      w = lattice.GetWannier(0, i);
      dw = lattice.GetLaplacian(0, i);
      ff[i] = -3.80998208024*w*dw + w*w*v;
    }
    return 0;
  };

  engine.divonne_result(integrand, b3, resN, errN, pN, steps);
  
  for (size_t i  = 0; i < resN.size(); ++i)
    std::cout << resN[i] << " "
              << errN[i] << " "
              << steps << "\n";




  return 0;

  /*
  std::ofstream file("orbital2.dat", std::ofstream::trunc);
  
  size_t npt = 110;
  for (size_t i  = 1; i < npt; ++i) {
    for (size_t j  = 1; j < npt; ++j) {
      for (size_t k  = 1; k < npt; ++k) {
        double x = -3.012 + 14.0/npt*i;
        double y = -6.012 + 11.0/npt*j;
        double z = -6.012 + 11.0/npt*k;

        lattice.UpdateWanniers({{x, y, z}});

        file << lattice.GetWannier(0, 12) << " ";
      }
    }
    std::cout << i << "\n";
  }
  
  file.close();
  */

  return 0;
}
