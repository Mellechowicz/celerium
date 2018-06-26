#include "../../headers/atomic_wave_function.h"
#include "../../headers/wfsolver_interface.h"
#include <iostream>
#include <iomanip>      // std::setprecision

using namespace celerium;

#define bohr_radius 0.529177249

// Manually define doubly differentiable radial wave function
// corresponding to the hydrogen 2p state.
// First and second derivatives are manually set
// through methods D1 and D2, respectively.
// This is not how one is normally intended to do that, but the object
// will serve us for testing purposes.
class RadialWaveFunction2p
{
   public:

  // Radial potential.
  double operator()(double r) {
    double c1 = 1.0/pow(2*bohr_radius, 3.0/2.0)/sqrt(3.0)/bohr_radius;
    double c2 = 0.5/bohr_radius;    
    return c1*r*exp(-c2*r);
  }

  // First derivative.
  double D1(double r) {
    double c1 = 1.0/pow(2.0*bohr_radius, 3.0/2.0)/sqrt(3.0)/bohr_radius;
    double c2 = 0.5/bohr_radius; 
    return c1*exp(-c2*r) - c1*c2*r*exp(-c2*r);
  }

  // Second derivative.
  double D2(double r) {
    double c1 = 1.0/pow(2.0*bohr_radius, 3.0/2.0)/sqrt(3.0)/bohr_radius;
    double c2 = 0.5/bohr_radius; 
    return -2*c1*c2*exp(-c2*r) + c1*c2*c2*r*exp(-c2*r);
  }
};

// Define discretized Coulomb potential for hydrogen atom.
// In realistic situation, some pseudopotential should be taken.
  class CoulombPotential{
    
   public:

    // Default configuration.
    CoulombPotential() {
      this->number_of_samples = 200;
      this->r_min = 0.0001;
      this->r_max = 30;
    }

    // Configure the Coulomb potential sampling.
    CoulombPotential(size_t number_of_samples, double r_min, double r_max) {
      this->number_of_samples = number_of_samples;
      this->r_min = r_min;
      this->r_max = r_max;
    }

    int get_mesh(std::vector<double> &v) const {
      v.clear();
      double x_min = log(r_min);
      double x_max = log(r_max);
      for (size_t i = 1; i < number_of_samples; ++i) {
        const double r =
            exp(x_min + (x_max - x_min)*i/((double)number_of_samples - 1.0));
        v.push_back(r);
      }
      return 0;
    }

    int get_local(std::vector<double> &v) const {
      v.clear();
      double x_min = log(r_min);
      double x_max = log(r_max);
      for (size_t i = 1; i < number_of_samples; ++i) {
        const double r =
            exp(x_min + (x_max - x_min)*i/((double)number_of_samples - 1.0));
        v.push_back(-14.399645352/r);
      }
      return 0;
    }

   private:
    size_t number_of_samples;
    double r_min;
    double r_max;
  };

int main()
{
  std::cout << std::fixed << std::setprecision(6);

  // Create radial potential object for 2p orbital.
  // It is manually defined and exact.
  RadialWaveFunction2p exact_radial_2p_wave_function;

  // Create complete wave function of the hydrogen 2p_z orbital (l=1, m=0).
  // At this point, it encompasses both radial and orbital parts.
  AtomicWaveFunction
      exact_2p_wave_function(exact_radial_2p_wave_function, 1, 0);

  // That's it. We can now evaluate the wave function and its laplacian:
  ArithmeticVector coords({0.1, 0.2, 0.5});
  std::cout << "\n\nValue of the manually defined 2p_z hydrogen wave function \
 and its laplacian for fixed coordinates:\n";
  std::cout << "Psi(0.1, 0.2, 0.5): "
            << exact_2p_wave_function(coords) << "\n";
  std::cout << "Laplacian Psi(0.1, 0.2, 0.5): "
            << exact_2p_wave_function.Laplacian(coords) << "\n";

  // In a realistic situation, the derivatives of the radial function are
  // are not explicitly known, e.g., for those frunctions obtained by 
  // numerically solving the Schroedinger equation.
  // We will now demonstrate how this case is handled. 
  // First, we need to determine numerically the doubly
  // differentiable radial function for the 2p orbital.

  // Create solver intended for the Coulomb potential object.
  WFSolver<CoulombPotential> solver;

  // Create an instance of CoulombPotential.
  CoulombPotential coulomb_potential;

  // Set the parameters for the solver.
  wf_solver_params params =
      {
        1e-3,
        0.01,
        1000,
        500,
        250,
        1e-8
      };

  // Now we obtain doubly differentiable radial wave function
  // for the Coulomb potential. The result will be stored
  // in the Interpolator object. Note that Inteprolator
  // is a doubly differentiable function. Everything will
  // be computed automatically and Laplacian method will
  // be available.
  double energy;
  Interpolator numerical_2p_radial_wave_function;
  int status = solver.GetEigenstate(2, 1,
                                    coulomb_potential,
                                    params,
                                    numerical_2p_radial_wave_function,
                                    energy);

  if (status == 0) {
    // Let us verify that the energy of the state is correct.
    std::cout << "\n\nLet us check if the energy of numerical solution for\
 2p orbital is correct:\n";
    
    std::cout << "abs(eigenenergy)/Ry*2^2: " << (fabs(energy)/13.605693009*4)
              << " (exact: 1.0)\n";
  }
  else {
    // Error, solver failed.
    throw std::runtime_error("WFSolver failed.");
  }

  // Create atomic wave function using calculated radial part.
  // Once again, l=1, m=0. This will be 2p_z orbital, cf.
  // https://en.wikipedia.org/wiki/Table_of_spherical_harmonics.
  AtomicWaveFunction<Interpolator>
      numerical_2p_wave_function(numerical_2p_radial_wave_function, 1, 0);

  // Compare the exact wave function with the numerically calculated one
  // and print the results.

  ArithmeticVector r0({1e-2, 1e-3, 1e-2});
  ArithmeticVector dr({0.4, 0.1, 0.3});

  std::cout << "\n\n";
  std::cout << "Comparison of the wave function values for the 2p_z prbital\n";
  std::cout << "along selected cut in space. Numerical solution\n";
  std::cout << "of the Schroedinger equation is";
  std::cout << " compared with analytical result.\n\n";
  
  std::cout << std::setw(9) << std::left << "x"
            << std::setw(9) << std::left << "y"
            << std::setw(15) << std::left << "z"
            << std::setw(15) << std::left << "Psi (exact)"
            <<  std::setw(15) << std::left << "Psi (numerical)\n";
  
  for (size_t i = 0; i < 15; ++i) {
    std::cout << std::setw(9) << std::left << (r0 + dr*i)[0]
              << std::setw(9) << std::left << (r0 + dr*i)[1]
              << std::setw(15) << std::left << (r0 + dr*i)[2]
              << std::setw(15) << std::left << exact_2p_wave_function(r0 + dr*i)
              << std::setw(15) << std::left << numerical_2p_wave_function(r0 + dr*i) << "\n";
  }
  
  std::cout << "\n\n";
    std::cout << "\n\n";
  std::cout << "Comparison of the wave function laplacian for the 2p_z prbital";
  std::cout << "\nalong selected cut in space. Numerical solution\n";
  std::cout << "of the Schroedinger equation is";
  std::cout << " compared with analytical result.\n\n";
  
  std::cout << std::setw(9) << std::left << "x" 
            << std::setw(9) << std::left << "y"
            << std::setw(21) << std::left << "z"
            << std::setw(21) << std::left << "Delta Psi (exact)"
            << std::setw(21) << std::left << "Delta Psi (numerical)\n";

  for (size_t i = 0; i < 15; ++i) {

    std::cout << std::setw(9) << std::left << (r0 + dr*i)[0]
              << std::setw(9) << std::left << (r0 + dr*i)[1]
              << std::setw(21) << std::left << (r0 + dr*i)[2]
              << std::setw(21) << std::left << exact_2p_wave_function.Laplacian(r0 + dr*i)
              << std::setw(21) << std::left << numerical_2p_wave_function.Laplacian(r0 + dr*i) << "\n";
  }

  std::cout << "\n\nTEST FINISHED...\n\n";
  
  return 0;
}
