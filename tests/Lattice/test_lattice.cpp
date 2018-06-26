#include <lattice.h>
#include <iomanip>      // std::setprecision
#include <functional>

#define bohr_radius 0.529177249

using namespace celerium;

int main(int argc, char *argv[])
{

  // Manually define Coulomb potential.
  std::function coulomb_potential = [](double r) -> double {return -14.399645352/r;};

  // Manually define radial wave function of the Hydrogen 2p oritals.
  std::function radial_wf_2p = [](double r) {
    double c1 = 1.0/pow(2*bohr_radius, 3.0/2.0)/sqrt(3.0)/bohr_radius;
    double c2 = 0.5/bohr_radius;    
    return c1*r*exp(-c2*r);
  };
  
  OrbitalClass hydrogen_2p(radial_wf_2p, -13.605693009/2/2, 1);

  std::vector orbitals ({hydrogen_2p});
  
  Element hydrogen("H", coulomb_potential, orbitals);

  ElementaryCell<decltype(hydrogen)> lattice;


  lattice.AddSite(hydrogen, {{0.13, 0.42, 0.3}});
  lattice.AddSite(hydrogen, {{0.13, 0.42, 0.35}});

  std::cout << "\n\n" << lattice.NOrbitals() << "\n\n";

  std::vector<double> result;
  //lattice.EvaluateOrbitals({{0.5, 0.5, 0.5}}, result);

  std::cout << "\n\n";
  for (auto r : result) std::cout << r << " ";
  std::cout << "\n\n";
  
  return 0;
}
