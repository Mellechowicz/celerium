#include <orbital_class.h>
#include <functional>
#include <iomanip>

using namespace celerium;

#define bohr_radius 0.529177249

int main(int argc, char *argv[])
{

  // Manually define radial wave function of the Hydrogen 2p oritals.
  std::function radial_wf_2p = [](double r) -> double {
    double c1 = 1.0/pow(2*bohr_radius, 3.0/2.0)/sqrt(3.0)/bohr_radius;
    double c2 = 0.5/bohr_radius;    
    return c1*r*exp(-c2*r);
  };

  OrbitalClass hydrogen_2p(radial_wf_2p,       // Wave function.
                           -13.605693009/2/2,  // energy = - Ry/4
                           2,                  // n = 2
                           1);                 // l = 1

  // hydrogen_2p represents entire set of hydrogen 2p orbitals.
  // We can demonstrate that by printing values of p_x, p_y, and p_z
  // orbitals along selected axes.

  std::cout << "\n\n";
  
  std::cout << "Hydrogen p_z orbital along z axis: \n\n";
  std::cout << std::setw(10) << std::left << "z";
  std::cout << std::setw(10) << std::left << "Psi(0, 0, z)";
  std::cout << "\n";
      
  for (size_t i = 0; i < 10; ++i) {
    std::cout << std::setw(10) << std::left << i/2.0;
    std::cout << std::setw(10) << std::left << hydrogen_2p.Eval({{0, 0, i/2.0}}, 0);
    std::cout << "\n";
  }

  std::cout << "\n\n";
  
  std::cout << "Hydrogen p_x orbital along x axis: \n\n";
  std::cout << std::setw(10) << std::left << "x";
  std::cout << std::setw(10) << std::left << "Psi(x, 0, 0)";
  std::cout << "\n";
  
  for (size_t i = 0; i < 10; ++i) {
    std::cout << std::setw(10) << std::left << i/2.0;
    std::cout << std::setw(10) << std::left << hydrogen_2p.Eval({{i/2.0, 0, 0}}, 1);
    std::cout << "\n";
  }


  std::cout << "\n\n";
  
  std::cout << "Hydrogen p_y orbital along y axis: \n\n";
  std::cout << std::setw(10) << std::left << "y";
  std::cout << std::setw(10) << std::left << "Psi(0, y, 0)";
  std::cout << "\n";


  for (size_t i = 0; i < 10; ++i) {
    std::cout << std::setw(10) << std::left << i/2.0;
    std::cout << std::setw(10) << std::left << hydrogen_2p.Eval({{0, i/2.0, 0}}, -1);
    std::cout << "\n";
  }

  std::cout << "\n\nActive m values:\n\n";
  
  for (auto m : hydrogen_2p.GetActiveMValues()) {
    std::cout << m << "\n";
  }

  std::cout << "\n\n TEST FINISHED...\n\n";
  
  return 0;
}
