#include <iomanip>      // std::setprecision
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
  Element chromium("Cr", potential_cr, {});

  // Add 3d-like orbitals. Note that mesh is not passed as an argument here
  // as it is provided by LocalPotential object. Also, solver parameters
  // are not specified (default values are used).
  chromium.AddOrbitalClass(3,         // n = 3
                           2);        // l = 2
                           
  // Add 4s-like orbitals.
  chromium.AddOrbitalClass(4,         // n = 4
                           0);        // l = 0

  ElementaryCell elementary_cell;

  std::array<ArithmeticVector, 3> basis;
  basis[0] = {{2.91, 0.0, 0.0}};
  basis[1] = {{0.0, 2.91, 0.0}};
  basis[2] = {{0.0, 0.0, 2.91}};
  elementary_cell.SetBasis(basis);
  
  elementary_cell.AddSite("Cr(1)", chromium, {{0, 0, 0}});
  elementary_cell.AddSite("Cr(2)", chromium, {{2.91*0.5, 2.91*0.5, 2.91*0.5}});


  /////////////////////////////////////////////////////
  /////////////////////////////////////////////////////
  /////////////////////////////////////////////////////
  /////////////////////////////////////////////////////


      

  std::cout << "NObitals: " <<  elementary_cell.NOrbitals() << "\n\n";

  std::cout << "Normalization test:\n\n";
  
  std::function<int(const double *, double *)> integrand;
  integrand = [&](const double *xx, double *ff) {
    std::vector<double> result (12);
    elementary_cell.EvaluateOrbitals({{xx[0], xx[1], xx[2]}}, result);
    for (size_t i = 0; i < 12; ++i)
      ff[i] = result[i]*result[i];
    return 0;
  };

      
  cuba::Cuba engine(1e7,1e5,1e-5);
  
  std::vector<std::pair<double,double>> b3(3,std::make_pair(-10,10));
  std::vector<double> resN (12), errN (12), pN (12);  
  int steps = 0;

  engine.divonne_result(integrand,b3,resN,errN,pN,steps);

  for (size_t i  = 0; i < resN.size(); ++i)
    std::cout << "\nint: "
              << resN[i] << " "
              << errN[i] << " "
              << steps << "";

  std::cout << "\n\n\KONIEC...\n\n\n";

      
      
  return 0;
}

