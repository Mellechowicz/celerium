#include <iomanip>      // std::setprecision
#include <elementary_cell.h>
#include <potential.h>
#include <newcubawrapper.h>

using namespace celerium;

int main(int argc, char *argv[])
{

  cubacores(1, 1);

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

  


  std::cout << "\n\n" << elementary_cell.NOrbitals() << "\n\n";
  
  std::vector<double> result;

  //ArithmeticVector v({0.1, 0.2, 0.4});
  //elementary_cell.EvaluateOrbitals(v, result);

  //for (auto r : result) {
    //std::cout << r << "\n";
  //}

  std::cout << "\n\n";
  
  //elementary_cell.EvaluateLaplacians(v, result);

  //for (auto r : result) {
    //std::cout << r << "\n";
  //}

  /*
  size_t n_pt = 30;
  double a = 2.91;
  for (size_t i = 0; i < n_pt; ++i) {
    for (size_t j = 0; j < n_pt; ++j) {
      std::cout << a*i/n_pt << " ";
      std::cout << a*j/n_pt << " ";
      std::cout << elementary_cell.CellPotential({{a*i/n_pt, a*j/n_pt, 2.79/4}}, 10.0);
      std::cout << "\n";
    }
  }
  */


  

//elementary_cell.EvaluatePotentials(v, result);


/*
  std::cout << "\n\n";
  
  for (size_t i  = 0; i < elementary_cell.NSites(); ++i) {
    std::cout << result[i] << "\n";
  }
*/


    //std::vector<double> res;
    //elementary_cell.EvaluateOrbitals({{0.2, 0.1, 0.15}}, res);

  



  std::cout << chromium.GetRadialPotential()(0.1) << "\n";

  Element aa = chromium;

    std::cout << aa.GetRadialPotential()(0.1) << "\n";
             
  
    std::function<int(const double*, double*)> integrand = [&](const double xx[], double ff[]) -> int {
    std::vector<double> res (12);
    //elementary_cell.EvaluateOrbitals({{0.1, 0.1, 0.1}}, res);

    
    
    for(size_t j = 0; j < 1; ++j)  ff[j] = 0.1;
    return 0;
  };

  double ff[10];
  double xx[10];
      
  integrand(xx, ff);

  std::cout << "i : " << ff[0] << "\n";
  
  //return 0;
   
  cuba::Cuba engine(30000,1000,1e-3);

  
  std::vector<std::pair<double,double>> b3(3,std::make_pair(0,1));

  std::array<double,1> resN, errN, pN;
  int steps = 0;

  engine.suave_result(integrand,b3,resN,errN,pN,steps);


  std::cout << "aaaaaaaa";
  
  for (size_t i  = 0; i < resN.size(); ++i)
  std::cout << "\nint: "
            << resN[i] << " "
            << errN[i] << " "
            << pN[i] << " "
            << steps << "";

  std::cout << "\n\n";
  
  return 0;
}

