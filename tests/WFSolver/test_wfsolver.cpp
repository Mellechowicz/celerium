#include <wfsolverinterface.h>
#include <iostream>
#include <iomanip>      // std::setprecision

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

    int get_mesh(std::vector<double> &v) {
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

    int get_local(std::vector<double> &v) {
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



  class ArbitraryPotential{
    
   public:

    int get_mesh(std::vector<double> &v) {
      v =
      {
        0.0,
        1.0,
        2.0,
        3.0,
        4.0,
        5.0,
        6.0,
        7.0,
        8.0,
        9.0,
        10.0,
        15.0,
        20.0,
        30.0,
        40.0,
        100.0
      };
      return 0;
    }

    int get_local(std::vector<double> &v) {
      v =
      {
        -3.0,
        -5.0,
        -7.0,
        -8.0,
        -6.0,
        -4.0,
        -3.0,
        -2.0,
        -1.0,
        -0.5,
        -0.1,
        -0.01,
        -0.0001,
        -0.00001,
        -0.0000001,
        -0.0000000001
      };
      return 0;
    }
  };

using namespace celerium;

int main()
{
  std::cout << std::fixed << std::setprecision(12);

  CoulombPotential coulomb_potential_fine(10000, 1e-5, 40);
  CoulombPotential coulomb_potential_coarse(50, 1e-5, 40);
  ArbitraryPotential arbitrary_potential;

  wf_solver_params params;
  params.r_min = 1e-5;
  params.energy_step = 0.001;
  params.matrix_dim = 2000;
  params.grid_size = 5000;
  params.matching_index = 2500;
  params.energy_accuracy = 1e-10;

  WFSolver<CoulombPotential>
      solver_fine(coulomb_potential_fine, params);
  WFSolver<CoulombPotential>
      solver_coarse(coulomb_potential_coarse, params);
  WFSolver<ArbitraryPotential>
      solver_arbitrary(arbitrary_potential, params);
   
  std::function<double(double)>
      coulomb_potential_exact = [](double r) {return -14.399645352/r;};


  
  
  std::vector<double> wave_function;
  double energy;
  std::vector<double> mesh;

  std::cout << "Test #1: 1s state for hydrogen "
            << "atom (exact Coulomb potential).\n";
  std::cout << "         Using low-level solver.\n";

  eigenstate_struct eigenstate;
  int test1_status = SolveRadialSchroedingerEqn(coulomb_potential_exact,
                                      0,
                                      1e-12,
                                      40,
                                      -13.6,
                                      0.1,
                                      100000,
                                      50000,
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
  std::cout << "Test #2: 3d state for hydrogen "
            << "atom (fine Coulomb potential discretization).\n";
  std::cout << "         Using high-level solver.\n";
  
  int test2_status = solver_fine.GetEigenstate(3, 2, wave_function, energy);

  std::cout << solver_fine.GetLog() << "\n";

  if (test2_status == 0) {
    std::cout << "\nCalculated energy: " << energy << " (";
    std::cout << "exact: " << -1.51174366767 << ")\n";
    std::cout << "sqrt(abs(Ry/energy)): "
              << sqrt(fabs(13.605693009/energy)) << " (exact: 3)";
    std::cout << "\nTEST #2 PASSED!\n";
  }
  else {
        std::cout << "\nTEST #2 FAILED!\n";
  }

  std::cout << "\n\nTest #3: 3d state for hydrogen "
            << "atom (coarse Coulomb potential discretization).\n";
  std::cout << "         Using high-level solver.\n";
  
  int test3_status = solver_coarse.GetEigenstate(3, 2, wave_function, energy);

  std::cout << solver_coarse.GetLog() << "\n";

  if (test3_status == 0) {
    std::cout << "\nResult:\n";
    std::cout << "Calculated energy: " << energy << " (";
    std::cout << "exact: " << -1.51174366767 << ")\n";
    std::cout << "sqrt(abs(Ry/energy)): "
              << sqrt(fabs(13.605693009/energy)) << " (exact: 3)";

    coulomb_potential_coarse.get_mesh(mesh);
    std::cout << "\n\nWave function (format: r, \Psi(r)):\n";
    for (size_t i = 0; i < wave_function.size(); ++i) {
      std::cout << mesh[i] << " " << wave_function[i] << "\n";
    }

    std::cout << "\nTEST #3 PASSED!\n";
  }
  else {
        std::cout << "\nTEST #3 FAILED!\n";
  }
  
  
  std::cout << "\n\nTest #3: 3p-like state for "
            << "custom potential.\n";
  std::cout << "         Using high-level solver.\n";
  
  int test4_status =
      solver_arbitrary.GetEigenstate(3, 1, wave_function, energy);

  std::cout << solver_arbitrary.GetLog() << "\n";

  if (test4_status == 0) {
    std::cout << "\nResult:\n";
    std::cout << "Calculated energy: " << energy << "\n";
    std::cout << "sqrt(abs(Ry/energy)): "
              << sqrt(fabs(13.605693009/energy));

    
    arbitrary_potential.get_mesh(mesh);
    std::cout << "\n\nWave function (format: r, \Psi(r)):\n";
    for (size_t i = 0; i < wave_function.size(); ++i) {
      std::cout << mesh[i] << " " << wave_function[i] << "\n";
    }

    std::cout << "\nTEST #4 PASSED!\n";
  }
  else {
    std::cout << "\nTEST #4 FAILED!\n";
  }
  
  std::cout << "\n\nTEST FINISHED...\n\n";
  
  return 0;
}
