#include <lattice.h>
#include <potential.h>

using namespace celerium;

int main(int argc, char *argv[])
{  
  // Load pseudopotential.
  LocalPotential potential_o("O.UPF");
  LocalPotential potential_cu("Cu.UPF");
  LocalPotential potential_la("La.UPF");

  Element oxygen("O",            // Element name.
                 potential_o,    // Potential.
                 {});            // Do not add any orbital classes manually.

  Element copper("Cu",            // Element name.
                 potential_cu,    // Potential.
                 {});             // Do not add any orbital classes manually.

  Element lanthanum("La",            // Element name.
                    potential_la,    // Potential.
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
  basis[2] = {{0.0, 0.0, c}};

  ElementaryCell elementary_cell(basis);
  
  elementary_cell.AddSite("Cu", copper, {{0, 0, 0}});
  elementary_cell.AddSite("O(1)", oxygen, {{a*0.5, 0.0, 0.0}});
  elementary_cell.AddSite("O(2)", oxygen, {{0.0, a*0.5, 0.0}});
  elementary_cell.AddSite("O(3)", oxygen, {{0.0, 0.0, c*0.18623}});
  elementary_cell.AddSite("O(4)", oxygen, {{0.0, 0.0, -c*0.18623}});

  //elementary_cell.AddSite("La(1)", lanthanum, {{a*0.5, a*0.5, c*0.13911}});
  //elementary_cell.AddSite("La(2)", lanthanum, {{a*0.5, a*0.5, -c*0.13911}}); 

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


  for (size_t i  = 0; i < elementary_cell.NOrbitals(); ++i) {
    std::cout << elementary_cell.GetOrbitalDescription(i) << "\n";
  }
   
  lattice.LoadWannierDataFromFile("wanniers.dat");

  cuba::Cuba engine(1e6,1e5,1e-7);
  engine.parameters.nnew = 1e3;
  engine.parameters.flatness = 400;
  
  std::vector<std::pair<double,double>> b3(3,std::make_pair(-10,10));
  std::vector<double> resN (5), errN (5), pN (5);  
  int steps = 0;

  std::function<int(const double *, double *)> integrand;
  integrand = [&](const double *xx, double *ff) {
                
    double w, dw, v;

    //ArithmeticVector coords({xx[0], xx[1], xx[2]});
    //lattice.UpdateWanniers(coords);
    //lattice.UpdateLaplacians(coords);
    
    for (size_t i = 4; i <= 8; ++i) {
      w = lattice.EvaluateWannier(0, i, xx);
      dw = lattice.EvaluateLaplacian(0, i, xx);
      //w = lattice.GetWannier(0, 8);
      //dw = lattice.GetLaplacian(0, 8);
      v = lattice.EvaluateCrystalPotential(xx);
      //std::cout << v << "\n";
      ff[i-4] = -3.80998208024*w*dw + w*w*v;
    }
    return 0;
  };

  
  engine.divonne_result(integrand, b3, resN, errN, pN, steps);
  
  for (size_t i  = 0; i < resN.size(); ++i)
    std::cout << resN[i] << " "
              << errN[i] << " "
              << steps << "\n";
  



  return 0;


  
  std::ofstream file("orbital2.dat", std::ofstream::trunc);
  
  size_t npt = 110;
  for (size_t i  = 1; i < npt; ++i) {
    for (size_t j  = 1; j < npt; ++j) {
      for (size_t k  = 1; k < npt*0+2; ++k) {
        double x = -a + 2*a/npt*i;
        double y = -a + 2*a/npt*j;
        y = 0;
        double z = -c + 2*c/npt*j;

        //lattice.UpdateWanniers({{x, y, z}});

        double xx [] =  {x, y, z};
        file << x << " " <<  z << " " << lattice.EvaluateCrystalPotential(xx) << "\n";
      }
    }
    std::cout << i << "\n";
  }
  
  file.close();
  

  return 0;
}
