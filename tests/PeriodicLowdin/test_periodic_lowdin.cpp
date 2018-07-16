#include <periodic_lowdin.h>
#include <gslmatrixcomplex.h>
#include <iomanip>
#include <algorithm>

using namespace celerium;

int main(int argc, char *argv[])
{
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
      { {0, 1, 0}, gsl::Matrix(n_orbitals, {0.12, 0.08, 0.08, 0.12}) });

  orbital_overlaps.push_back(
      { {0, 2, 0}, gsl::Matrix(n_orbitals, {0.0, 0.0, 0.0, 0.0}) });

  orbital_overlaps.push_back(
      { {1, -2, 0}, gsl::Matrix(n_orbitals, {0.0, 0.0, 0.0, 0.0}) });

  orbital_overlaps.push_back(
      { {1, -1, 0}, gsl::Matrix(n_orbitals, {0.0, 0.0, 0.0, 0.0}) });

  orbital_overlaps.push_back(
      { {1, 0, 0}, gsl::Matrix(n_orbitals, {0.12, 0.08, 0.08, 0.12}) });

  orbital_overlaps.push_back(
      { {1, 1, 0}, gsl::Matrix(n_orbitals, {0.0, 0.3, 0.3, 0.0}) });

  orbital_overlaps.push_back(
      { {1, 2, 0}, gsl::Matrix(n_orbitals, {0.0, 0.0, 0.0, 0.0}) });

  orbital_overlaps.push_back(
      { {2, -2, 0}, gsl::Matrix(n_orbitals, {0.0, 0.0, 0.0, 0.0}) });

  orbital_overlaps.push_back(
      { {2, -1, 0}, gsl::Matrix(n_orbitals, {0.0, 0.0, 0.0, 0.0}) });

  orbital_overlaps.push_back(
      { {2, 0, 0}, gsl::Matrix(n_orbitals, {0.2, 0.0, 0.0, 0.2}) });

  orbital_overlaps.push_back(
      { {2, 1, 0}, gsl::Matrix(n_orbitals, {0.0, 0.0, 0.0, 0.0}) });

  orbital_overlaps.push_back(
      { {2, 2, 0}, gsl::Matrix(n_orbitals, {0.0, 0.0, 0.0, 0.0}) });


  PeriodicOthogonalization(orbital_overlaps,
                           {10, 10, 0},
                           {{0, 0, 0}, {2, 0, 0}},
                           orbital_positions,
                           wannier_coefficients);


  std::cout << "solution: " << "\n";
  
  for (size_t wannier_index = 0;
       wannier_index < wannier_coefficients[0].size();
       ++wannier_index) {
    
    std::cout << "Wannier #" << wannier_index << ":\n";

    for (size_t orbital_position_index = 0;
         orbital_position_index < orbital_positions.size();
         ++orbital_position_index) {

      std::cout << "(";
      std::cout << orbital_positions[orbital_position_index][0] << " ";
      std::cout << orbital_positions[orbital_position_index][1] << " ";
      std::cout << orbital_positions[orbital_position_index][2] << "): ";
      
      for (size_t orbital_index = 0;
           orbital_index < n_orbitals;
           ++orbital_index) {
        std::cout <<
            wannier_coefficients[0]
                                [wannier_index]
                                [orbital_index+n_orbitals*orbital_position_index] << " ";
      }
      std::cout << "\n";
    }
    std::cout << "\n";
  }

  std::cout << "\n\n";
  
  std::cout << "Wannier orbital_overlaps: \n";

  for (size_t o1 = 0; o1 < 2; ++o1) {
    for (size_t o2 = 0; o2 < 2; ++o2) {
      std::cout << ScalarProduct(orbital_overlaps,
                                 orbital_positions,
                                 wannier_coefficients,
                                 false,
                                 0, o1,
                                 0, o2) << " ";
    }
    std::cout << "\n";
  }
  
  

      


  return 0;
}
