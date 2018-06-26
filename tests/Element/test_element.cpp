#include <cmath>
#include <functional>
#include <element.h>
#include <orbital_class.h>
#include <potential.h>
#include <wfsolver_interface.h>

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

  // Manually define radial wave function of the hydrogen 2p oritals.
  std::function radial_wf_2p = [](double r) {
    double c1 = 1.0/pow(2*bohr_radius, 3.0/2.0)/sqrt(3.0)/bohr_radius;
    double c2 = 0.5/bohr_radius;    
    return c1*r*exp(-c2*r);
  };
  
  // Create the class of hydrigen 2p orbitals (by default it
  // encompasses states with m = -1, 0, +1).
  OrbitalClass hydrogen_2p(radial_wf_2p,        // wave function 
                           -13.605693009/2/2,   // energy
                           1);                  // l = 1

  // Create hydrogen with defined 2p orbitals.
  Element hydrogen("H", coulomb_potential, {hydrogen_2p});

  // *****************************************************
  // * USE CASE 2: manually define analytical potential. *
  // * Wave functions are calclated numerically.         * 
  // *****************************************************

  // ***********************************************
  // * USE CASE 3: load pseudopotential from file. *
  // * Wave functions are calclated numerically.   *
  // ***********************************************

  LocalPotential potential_cr("../Potential/Cr.UPF");

  // Initialize element using chromium pseudopotential.
  Element chromium("Cr", potential_cr);

  // Set up the solver.
  wf_solver_params params;
  params.r_min = 1e-6;
  params.energy_step = 0.01;
  params.matrix_dim = 2000;
  params.grid_size = 10000;
  params.matching_index = 5000;
  params.energy_accuracy = 1e-8;

  /*
  // Add 3d-like orbitals.
  chromium.AddOrbitalClass(3,         // n = 3
                           2,         // l = 2
                           params);   // solver parameters

  // Add 4s-like orbitals.
  chromium.AddOrbitalClass(4,         // n = 4
                           0,         // l = 0
                           params);   // solver parameters
  
  
  for (auto &orbital_class : chromium.GetOrbitalClasses()) {
    std::cout << orbital_class.Eval({{0.1, 0.2, 0.3}}, 0) << "\n";
  }
  */

  
 
  return 0;
}


