#include <iomanip>           // std::setprecision
#include <elementary_cell.h>
#include <potential.h>
#include <functional>
#include <newcubawrapper.h>

using namespace celerium;

int main(int argc, char *argv[])
{
  
  // Load pseudopotential.
  LocalPotential potential_cr("../Potential/Cr.UPF");

  // Initialize element using chromium pseudopotential.
  Element chromium("Cr",            // Element name.
                   potential_cr,    // Potential.
                   {});             // Do not add any orbital classes manually.

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

  /////////////////////////////////////////////////////

  std::cout << "Number of orbitals in the elementary cell: "
            <<  elementary_cell.NOrbitals() << "\n\n";

  std::vector<std::string> orbital_descriptions;
  elementary_cell.GetOrbitalDescriptions(orbital_descriptions);

  std::cout << "Description of the orbitals: \n";
  for (const auto &orbital_description : orbital_descriptions) {
    std::cout << orbital_description << "\n";
  }


  elementary_cell.SetCrystalPotentialCutoff(40.0);

    size_t npt = 100;
   for (size_t i  = 1; i < npt; ++i) {
     for (size_t j  = 1; j < npt; ++j) {
      
            double x = -5.012 + 10.024/npt*i;
            double y = -5.012 + 10.024/npt*j;
            double r [] = {x, y, a/4};
            std::cout << x << " " << y << " " <<
              elementary_cell.EvaluateCrystalPotential(r) << "\n";
    }
  }  

  return 0;

  std::cout << "\nOrbital normalization test:\n";
  
  cuba::Cuba engine(1e8,1e5,1e-4);
  
  std::vector<std::pair<double,double>> b3(3,std::make_pair(-10,10));
  std::vector<double> resN (12), errN (12), pN (12);  
  int steps = 0;

  std::function<int(const double *, double *)> integrand;
  integrand = [&](const double *xx, double *ff) {
    elementary_cell.EvaluateOrbitals(xx, ff);
    for (size_t i = 0; i < 12; ++i) ff[i] *= ff[i];
    return 0;
  };

  engine.divonne_result(integrand, b3, resN, errN, pN, steps);
  
  for (size_t i  = 0; i < resN.size(); ++i)
    std::cout << resN[i] << " "
              << errN[i] << " "
              << steps << "\n";

  std::cout << "\n\nTest finished...\n\n";
  
  return 0;
}

