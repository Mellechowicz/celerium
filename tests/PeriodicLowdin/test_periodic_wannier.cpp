#include <periodic_wannier.h>
#include <gslmatrixcomplex.h>
#include <iomanip>
#include <algorithm>
#include <sstream>

using namespace celerium;

int main(int argc, char *argv[])
{
  std::cout << std::fixed;
  std::cout << std::setprecision(12);
  
  // This object will store the contributions of respective orbitals to
  // Wannier functions.
  std::vector<std::vector<std::vector<double>>> wannier_coefficients;

  // This vector will be filled with the positions of orbitals
  // contributing to the Wannier function.
  std::vector<std::array<int, 3>> orbital_positions;

  // Overlap integrals between orbitals in different lattice positions.
  // Those are provided by the user.
  std::vector<std::pair<std::array<int, 3>, gsl::Matrix>> orbital_overlaps;

  size_t n_orbitals = 2;
  
  orbital_overlaps.push_back(
      { {0, 0, 0}, gsl::Matrix(n_orbitals, {1.0, 0.1, 0.1, 1.0}) });

  orbital_overlaps.push_back(
      { {0, 1, 0}, gsl::Matrix(n_orbitals, {0.2, 0.02, 0.01, 0.2}) });

  orbital_overlaps.push_back(
      { {0, 2, 0}, gsl::Matrix(n_orbitals, {0.0, 0.0, 0.0, 0.0}) });

  orbital_overlaps.push_back(
      { {1, -2, 0}, gsl::Matrix(n_orbitals, {0.0, 0.0, 0.0, 0.0}) });

  orbital_overlaps.push_back(
      { {1, -1, 0}, gsl::Matrix(n_orbitals, {0.0, 0.0, 0.0, 0.0}) });

  orbital_overlaps.push_back(
      { {1, 0, 0}, gsl::Matrix(n_orbitals, {0.2, 0.02, 0.01, 0.2}) });

  orbital_overlaps.push_back(
      { {1, 1, 0}, gsl::Matrix(n_orbitals, {0.0, 0.0, 0.0, 0.0}) });

  orbital_overlaps.push_back(
      { {1, 2, 0}, gsl::Matrix(n_orbitals, {0.0, 0.0, 0.0, 0.0}) });

  orbital_overlaps.push_back(
      { {2, -2, 0}, gsl::Matrix(n_orbitals, {0.0, 0.0, 0.0, 0.0}) });

  orbital_overlaps.push_back(
      { {2, -1, 0}, gsl::Matrix(n_orbitals, {0.0, 0.0, 0.0, 0.0}) });

  orbital_overlaps.push_back(
      { {2, 0, 0}, gsl::Matrix(n_orbitals, {0.0, 0.0, 0.0, 0.0}) });

  orbital_overlaps.push_back(
      { {2, 1, 0}, gsl::Matrix(n_orbitals, {0., 0.0, 0.0, 0.0}) });

  orbital_overlaps.push_back(
      { {2, 2, 0}, gsl::Matrix(n_orbitals, {0.0, 0.0, 0.0, 0.0}) });


  WannierData wannier_data;

  std::array<ArithmeticVector, 3> basis;
  basis[0] = ArithmeticVector({1, 0, 0});
  basis[1] = ArithmeticVector({0, 1, 0});
  basis[2] = ArithmeticVector({0, 0, 1});
  
  wannier_data.Initialize(orbital_overlaps, {{3, 3, 0}}, {{0, 0, 0}}, basis, 0.1);

  for (const auto& orbital : wannier_data.GetOrbitalPositions()) {
    std::cout << orbital.absolute_position << "\n";
  }

  
  for (const auto& wannier : wannier_data.GetWanniers()) {
    std::cout << "Wannier #" << wannier.wannier_index << ":\n";
    for (const auto &extended_coeff : wannier.extended_coeffs) {
      std::cout << "(";
      std::cout <<
          wannier_data.GetOrbitalPositions()[extended_coeff.orbital_index].cell_position[0] << ", ";
      std::cout <<
          wannier_data.GetOrbitalPositions()[extended_coeff.orbital_index].cell_position[1] << ", ";
      std::cout <<
          wannier_data.GetOrbitalPositions()[extended_coeff.orbital_index].cell_position[2] << "), ";

      std::cout << extended_coeff.orbital_index << ": ";
      std::cout << extended_coeff.coeff << "\n";
    }
    std::cout << "\n\n";
  }


  std::cout << "\n\n";

  /*
  PeriodicOthogonalization(
      orbital_overlaps,
      {{5, 5, 0}},             // Range of orbitals contributing to Wannier functions.
      {{0, 0, 0}, {1, 0, 0}},  // Requested positions of Wannier functions.
      orbital_positions,
      wannier_coefficients);
  */




  /*
  std::cout << "\n\n";
  
  for (size_t wannier_index = 0;
       wannier_index < wannier_coefficients[0].size();
       ++wannier_index) {
    
    std::cout << "Wannier #" << wannier_index << ":\n";
    std::cout << std::setw(20) << std::left << "Orb. position";
    std::cout << std::setw(20) << std::left << "Orb. #0 weight";
    std::cout << std::setw(20) << std::left << "Orb. #1 weight";
    std::cout << "\n";

    for (size_t orbital_position_index = 0;
         orbital_position_index < orbital_positions.size();
         ++orbital_position_index) {

      std::stringstream ss;
      ss << "(";
      ss << orbital_positions[orbital_position_index][0] << " ";
      ss << orbital_positions[orbital_position_index][1] << " ";
      ss << orbital_positions[orbital_position_index][2] << "):";
      
      std::cout << std::setw(20) << std::left << ss.str();
      
      for (size_t orbital_index = 0;
           orbital_index < n_orbitals;
           ++orbital_index) {
        std::cout << std::setw(20) <<
            wannier_coefficients[0]
                                [wannier_index]
                                [orbital_index+n_orbitals*orbital_position_index] << " ";
      }
      std::cout << "\n";
    }
    std::cout << "\n";
  }

  std::cout << "\n\n";
  
  std::cout << "Wannier overlap matrix: \n\n";

  std::cout << "On-site:\n";
  
  for (size_t o1 = 0; o1 < 2; ++o1) {
    for (size_t o2 = 0; o2 < 2; ++o2) {
      std::cout << ScalarProduct(orbital_overlaps,
                                 orbital_positions,
                                 wannier_coefficients,
                                 0, o1,
                                 0, o2) << " ";
    }
    std::cout << "\n";
  }
  
  std::cout << "\nDifferent sites:\n";
  
  for (size_t o1 = 0; o1 < 2; ++o1) {
    for (size_t o2 = 0; o2 < 2; ++o2) {
      std::cout << ScalarProduct(orbital_overlaps,
                                 orbital_positions,
                                 wannier_coefficients,
                                 0, o1,
                                 1, o2) << " ";
    }
    std::cout << "\n";
  }
        
  std::cout << "\n\nTest finished...\n\n";
  */
  return 0;
}
