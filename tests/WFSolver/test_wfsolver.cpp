#include <wfsolver.h>
#include <iostream>
#include <iomanip>      // std::setprecision

#include <functional>

using namespace celerium;

class CoulombPotential{
    
 public:

  CoulombPotential() {
    this->number_of_samples = 100;
    this->r_min = 0.0001;
    this->r_max = 50;
  }
    
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


int main() {

  std::cout << std::setprecision(12);

  std::cout << "Test #1: 1s state for hydrogen "
            << "atom (exact Coulomb potential).\n";
  std::cout << "         Using low-level solver.\n";

  std::function potential = [](double r) {return -14.399645352/r;};
  eigenstate_struct eigenstate;

  int test1_status = SolveRadialSchroedingerEqn(potential,
                                                0,
                                                1e-8,
                                                40.0,
                                                -13.5,
                                                0.01,
                                                10000,
                                                5000,
                                                1e-10,
                                                eigenstate);
  if (test1_status == 0) {
    std::cout << "\nCalculated energy: " << eigenstate.energy << " (";
    std::cout << "exact: " << -13.605693009 << ")\n";
    std::cout << "sqrt(abs(Ry/energy)): "
              << sqrt(fabs(13.605693009/eigenstate.energy)) <<" (exact: 1)\n";
    std::cout << "\nTEST #1 PASSED!\n";
  }
  else {
    std::cout << "\nTEST #1 FAILED!\n";
  }

  std::cout << "\n\n";
  
  std::cout << "Test #2: 1s state for hydrogen "
            << "atom (coarsly inteporlated Coulomb potential).\n";
  std::cout << "         Using high-level solver.\n";

  WFSolver solver;
  CoulombPotential coulomb_potential_coarse(50, 1e-6, 40);
  std::vector<double> mesh;
  std::vector<double> values;
  coulomb_potential_coarse.get_mesh(mesh);
  coulomb_potential_coarse.get_local(values);
  
  int test2_status =   solver.GetEigenstate(1, 0,
                                            mesh,
                                            values,
                                            mesh,
                                            eigenstate);
  
  if (test2_status == 0) {
    std::cout << "\nCalculated energy: " << eigenstate.energy << " (";
    std::cout << "exact: " << -13.605693009 << ")\n";
    std::cout << "sqrt(abs(Ry/energy)): "
              << sqrt(fabs(13.605693009/eigenstate.energy)) <<" (exact: 1)\n";
    std::cout << "\nTEST #2 PASSED!\n";
  }
  else {
    std::cout << "\nTEST #2 FAILED!\n";
  }

  std::cout << "\n\n";
  
  std::cout << "Test #3: 1s state for hydrogen "
            << "atom (fine inteporlation of the Coulomb potential).\n";
  std::cout << "         Using high-level solver.\n";

  CoulombPotential coulomb_potential_fine(10000, 1e-6, 40);
  coulomb_potential_fine.get_mesh(mesh);
  coulomb_potential_fine.get_local(values);
  
  int test3_status =   solver.GetEigenstate(1, 0,
                                            mesh,
                                            values,
                                            mesh,
                                            eigenstate);
  
  if (test3_status == 0) {
    std::cout << "\nCalculated energy: " << eigenstate.energy << " (";
    std::cout << "exact: " << -13.605693009 << ")\n";
    std::cout << "sqrt(abs(Ry/energy)): "
              << sqrt(fabs(13.605693009/eigenstate.energy)) <<" (exact: 1)\n";
    std::cout << "\nTEST #3 PASSED!\n";
  }
  else {
    std::cout << "\nTEST #3 FAILED!\n";
  }

  std::cout << "\n\nTest finished...\n\n";

  return 0;
};


