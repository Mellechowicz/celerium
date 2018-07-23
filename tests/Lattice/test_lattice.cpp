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

 


  cuba::Cuba engine(5*1e5,5*1e5,1e-4);
  engine.parameters.flatness = 3000;
  engine.parameters.nnew = 3000;
  engine.parameters.nmin = 1000;


  std::vector<std::pair<double, double>> integration_limits(3, {-15, 15});
  

  /*
  lattice.CalculateWannierCoefficients({{2, 2, 2}},
                                       {{4, 4, 4}},
                                       {{0, 0, 0}, {1, 0, 0}},
                                       integration_limits,
                                       engine,
                                       0.0001,
                                       true);
  
  lattice.SaveWannierDataToFile("data.dat");

  return 0;
  */
  lattice.LoadWannierDataFromFile("data.dat");


  size_t nmax = 12;
  
  std::function<int(double*, double*)> overlap =
      [&](double *xx, double *ff) {
    double wanniers_r1 [12];
    double wanniers_r2 [12];
    
    lattice.UpdateWanniers({{xx[0], xx[1], xx[2]}});
    lattice.GetWanniers(0, wanniers_r1);
    
    lattice.UpdateWanniers({{xx[3], xx[4], xx[5]}});
    lattice.GetWanniers(0, wanniers_r2);
    
    for (size_t i = 0; i < nmax; ++i) {
      ff[i] = wanniers_r1[i]*wanniers_r1[i]*wanniers_r2[i]*wanniers_r2[i];
    double dr =
        (xx[0]-xx[3])*(xx[0]-xx[3]) + 
        (xx[1]-xx[4])*(xx[1]-xx[4]) +
        (xx[2]-xx[5])*(xx[2]-xx[5]);
    dr = std::sqrt(dr);
    ff[i] *= 14.399645352/dr;
    }
    return 0;
  };

  std::vector<double> resN (nmax);
  std::vector<double> errN (nmax);
  std::vector<double> pN (nmax);
  int steps;

  std::vector<std::pair<double,double>> b6(6,std::make_pair(-10,10));  
  engine.suave_result(overlap, b6, resN, errN, pN, steps);

    size_t l {0};
    for (size_t i = 0; i < nmax; ++i) {
      std::cout << resN[l] << "     " << errN[l] << "   " << steps << "\n";
      ++l;
    }


  std::cout << "\n";
  
  return 0;





  
  /*
  size_t nmax = 12;
  
  std::function<int(double*, double*)> overlap =
      [&](double *xx, double *ff) {
    lattice.UpdateWanniers({{xx[0], xx[1], xx[2]}});
    double wanniers [12];
    lattice.GetWanniers(0, wanniers);
    size_t l {0};
    for (size_t i = 0; i < nmax; ++i) {
      for (size_t j = 0; j < nmax; ++j) {
        ff[l] = wanniers[i]*wanniers[j];
        ++l;
      }
    }
    return 0;
  };

  std::vector<double> resN (nmax*nmax);
  std::vector<double> errN (nmax*nmax);
  std::vector<double> pN (nmax*nmax);
  int steps;
  std::vector<std::pair<double,double>> b3(3,std::make_pair(-10,10));
  
  engine.suave_result(overlap, b3, resN, errN, pN, steps);

    size_t l {0};
    for (size_t i = 0; i < nmax; ++i) {
      for (size_t j = 0; j < nmax; ++j) {
        std::cout << resN[l] << "     " << errN[l] << "   " << steps << "\n";
        ++l;
      }

      std::cout << "\n\n\n";
    }

  std::cout << "\n";
  
  return 0;
*/  

  
  /*
  size_t npt = 100;
  for (size_t i  = 1; i < npt; ++i) {
    for (size_t j  = 1; j < npt; ++j) {
          for (size_t k  = 1; k < npt; ++k) {
      

            double x = -5.012 + 10.0/npt*i;
            double y = -5.012 + 10.0/npt*j;
            double z = -5.012 + 10.0/npt*k;
      
            lattice.UpdateWanniers({{x, y, z}});
            // std::cout << x << " " << y << " " << lattice.GetWannier(0, 11) << "\n";
          }
    }
  }
  */
  return 0;
}
