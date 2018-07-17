#include <cmath>
#include <functional>
#include <iomanip>
#include <element.h>
#include <orbital_class.h>
#include <potential.h>


using namespace celerium;

#define bohr_radius 0.529177249

int main(int argc, char *argv[])
{
  
  // *****************************************************************
  // * USE CASE 1: manually define everything.                       *
  // * Analytical form of the potential and wave functions is known. *
  // *****************************************************************
   
  // Manually define Coulomb potential.
  std::function coulomb_potential = [](double r) {return -14.399645352/r;};

  // Manually define radial wave function for the hydrogen 2p oritals.
  std::function radial_wf_2p = [](double r) {
    double c1 = 1.0/pow(2*bohr_radius, 3.0/2.0)/sqrt(3.0)/bohr_radius;
    double c2 = 0.5/bohr_radius;    
    return c1*r*exp(-c2*r);
  };
  
  // Create the class of hydrogen 2p orbitals (by default it
  // encompasses states with m = -1, 0, +1).
  OrbitalClass hydrogen_2p(radial_wf_2p,        // wave function 
                           -13.605693009/2/2,   // energy
                           2,                   // n = 2
                           1);                  // l = 1

  // Create hydrogen with 2p orbitals.
  Element hydrogen_manual("H", coulomb_potential, {hydrogen_2p});

  std::cout << "Values of wave functions along the " <<
      "selected axis for hydrogen 2p orbitals (USE CASE 1: evertyhing manually defined).\n\n";
  std::cout << std::setprecision(5);
  for (auto &orbital_class : hydrogen_manual.GetOrbitalClasses()) {
    std::cout << std::setprecision(15);
    std::cout << "n = " << orbital_class.GetN()
              << ", l = " << orbital_class.GetL() << ", E = "
              << orbital_class.GetEnergy() << "\n";
    std::cout << std::setprecision(5);
    for (auto m : orbital_class.GetActiveMValues()) {
      std::cout << "m = " << m << ", ";
      for (int i = 0; i < 10; ++i) 
        std::cout << std::setw(9) << orbital_class.Eval({{0.1, 0.3, i*0.5}}, m);
      std::cout << "\n";
    }
  }

  std::cout << "\n\n";
  
  // *****************************************************
  // * USE CASE 2: potential is known analytically.      *
  // * Wave functions are calclated numerically.         * 
  // *****************************************************

  // We first create hydrogen atom without any orbitals.
  // Those will be determined numerically.
  Element hydrogen_numerical("H", coulomb_potential, {});

  // Define the mesh over which the wave function will be sampled.
  // We need to do that as the std::function does not provide any mesh.
  std::vector<double> result_mesh;
  for (int i  = 0; i < 100; ++i) result_mesh.push_back(i*0.3);

  hydrogen_numerical.AddOrbitalClass(2,               // n = 2
                                     1,               // l = 1
                                     result_mesh);    // wave function sampling mesh

  std::cout << "Values of wave functions along the " <<
      "selected axis for hydrogen 2p orbitals (USE CASE 2: manually defined " <<
      "potential, numerical wave function).\n\n";
  std::cout << std::setprecision(5);
  for (auto &orbital_class : hydrogen_numerical.GetOrbitalClasses()) {
    std::cout << std::setprecision(15);
    std::cout << "n = " << orbital_class.GetN()
              << ", l = " << orbital_class.GetL() << ", E = "
              << orbital_class.GetEnergy() << "\n";
    std::cout << std::setprecision(5);
    for (auto m : orbital_class.GetActiveMValues()) {
      std::cout << "m = " << m << ", ";
      for (int i = 0; i < 10; ++i) 
        std::cout << std::setw(9) << orbital_class.Eval({{0.1, 0.3, i*0.5}}, m);
      std::cout << "\n";
    }
  }
  
  std::cout << "\n\n";
  
  // ***********************************************
  // * USE CASE 3: load pseudopotential from file. *
  // * Wave functions are calclated numerically.   *
  // ***********************************************

  // Load pseudopotential.
  LocalPotential potential_cr("../Potential/Cr.UPF");

  // Initialize element using chromium pseudopotential.
  // Do not add any orbitals manually.
  Element chromium("Cr", potential_cr, {});
  
  // Add 3d-like orbitals. Note that mesh is not passed as an argument here
  // as it is provided by LocalPotential object.
  chromium.AddOrbitalClass(3,         // n = 3
                           2);        // l = 2


  // Add 4s-like orbitals.
  chromium.AddOrbitalClass(4,         // n = 4
                           0);        // l = 0
  
  
  std::cout << "Values of wave functions along the " <<
      "selected axis for chromium 3d and 4s orbitals (USE CASE 3: full numerics).\n\n";
  std::cout << std::setprecision(5);
  for (auto &orbital_class : chromium.GetOrbitalClasses()) {
    std::cout << std::setprecision(15);
    std::cout << "n = " << orbital_class.GetN()
              << ", l = " << orbital_class.GetL() << ", E = "
              << orbital_class.GetEnergy() << "\n";
    std::cout << std::setprecision(5);
    for (auto m : orbital_class.GetActiveMValues()) {
      std::cout << "m = " << m << ", ";
      for (int i = 0; i < 5; ++i) 
        std::cout << std::setw(15) << orbital_class.Eval({{0.1, 0.3, i*0.5}}, m);
      std::cout << "\n";
    }
  }
    
  std::cout << "\n\nTEST FINISHED...\n\n";
 
  return 0;
}


