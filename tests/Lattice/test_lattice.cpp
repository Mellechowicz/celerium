#include <lattice.h>
#include <potential.h>

using namespace celerium;

  int fun(double *xx, double *ff) {
    ff[0] = 0;
    ff[1] = 0;
    return 0;
  }



int main(int argc, char *argv[])
{
  /*
  cuba::Cuba integrator(1e2, 1e1, 1e-4);
  cuba::Cuba integrator2(1e2, 1e1, 1e-4);

  std::vector<std::pair<double,double>> b3(3,std::make_pair(-10,10));
  std::vector<std::pair<double,double>> b32(2,std::make_pair(-1,1));
  std::vector<double> resN (2);
  std::vector<double> errN (2);
  std::vector<double> pN (2);
  int steps;

  std::function<int(double*,double*)> integrand;
  std::function<int(double*,double*)> integrand2;
  
  integrand = [](double *xx, double *ff) {
    ff[0] = 1;
    ff[1] = 1;
    //std::cout << "a";
    return 0;
  };

  integrand2 = [](double *xx, double *ff) {
    ff[0] = 0;
    ff[1] = 0;
    //std::cout << "b";
    return 0;
  };

  integrator.suave_result(integrand, b3, resN, errN, pN, steps);
  std::cout << "int " << resN[0] << " " << resN[1] << "\n";
  resN[0] = 0; resN[1] = 0;
  integrator.divonne_result(integrand2, b32, resN, errN, pN, steps);
    
  std::cout << "int " << resN[0] << " " << resN[1] << "\n";

  return 0;
  
  //cubacores(1, 1);
  */
  
  // Load pseudopotential.
  LocalPotential potential_cr("../Potential/Cr.UPF");

  // Initialize element using chromium pseudopotential.
  Element chromium("Cr",            // Element name.
                   potential_cr,    // Potential.
                   {});             // Do not add any orbital classes manually.

  // Now add 3d-like orbitals. Note that mesh is not passed as an argument here
  // as it is provided by LocalPotential object. Also, solver parameters
  // are not specified (default values are used).
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

 
  //cuba::Cuba engine(1e9,1e4,1e-4);

  cuba::Cuba engine(1e7,1e5,1e-4);
  
  //lattice.CalculateWannierCoefficients({{2, 2, 2}}, {{2, 2, 2}}, {{0, 0, 0}}, engine, true);
  //lattice.SaveWannierDataToFile("data.dat");

  //return 0;

  lattice.LoadWannierDataFromFile("data.dat");

  //lattice.PrintWannier(5, 0);

  
  std::function<int(double*, double*)> overlap =
      [&](double *xx, double *ff) {
    lattice.UpdateWanniers({{xx[0], xx[1], xx[2]}});
    ff[0] = lattice.GetWannier(0, 0)*lattice.GetWannier(0, 0);
    return 0;
  };

  std::vector<double> resN (1);
  std::vector<double> errN (1);
  std::vector<double> pN (1);
  int steps;
  std::vector<std::pair<double,double>> b3(3,std::make_pair(-10,10));
  
  engine.divonne_result(overlap, b3, resN, errN, pN, steps);

  std::cout << resN[0] << " ";
  std::cout << "\n";
  
  return 0;
  

  

  size_t npt = 100;
  for (size_t i  = 1; i < npt; ++i) {
    for (size_t j  = 1; j < npt; ++j) {

      double x = -5.012 + 10.0/npt*i;
      double y = -5.012 + 10.0/npt*j;
      
      lattice.UpdateWanniers({{x, y, 1.45}});
      std::cout << x << " " << y << " " << lattice.GetWannier(5, 0) << "\n";

    /*
    std::vector<double> result;
    lattice.EvaluateOrbitals({{0.1*i, 0.00, 0.00}},
                     {{0, 0, 0}},
                     result);

    std::cout << result[5] << "\n";
    */
    }
  }
  return 0;
}
