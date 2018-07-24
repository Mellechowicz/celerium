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
  
  wannier_data.Initialize(orbital_overlaps,
                          {{10, 10, 0}},
                          {{0, 0, 0}, {1, 0, 0}},
                          basis,
                          0.00);
  
  for (const auto& wannier : wannier_data.GetWanniers()) {
    std::cout << "Wannier #" << wannier.wannier_index << ":\n";
      std::cout << std::setw(20) << std::left << "Cell coords";
      std::cout << std::setw(20) << std::left << "Orbital index";
      std::cout << std::setw(20) << std::left << "Weight";
      std::cout << "\n";
    
    for (const auto &extended_coeff : wannier.extended_coeffs) {

      const auto &orbital_basis_coords =
          wannier_data.GetOrbitalPositions()
              [extended_coeff.position_index].basis_coords;
      
      std::stringstream ss;
      ss << "(" << orbital_basis_coords[0] << ", " <<
                   orbital_basis_coords[1] << ", " <<
                   orbital_basis_coords[2] << ")";
      
      std::cout << std::setw(20) << std::left << ss.str();
      std::cout << std::setw(20) << std::left << extended_coeff.orbital_index;
      std::cout << std::setw(20) << std::left << extended_coeff.coeff;
      std::cout << "\n";
    }
    std::cout << "\n\n";
  }


  std::cout << "Wannier overlaps (on-site:)\n";

  for (size_t i  = 0; i < n_orbitals; ++i) {
    for (size_t j  = 0; j < n_orbitals; ++j) {
      std::cout << ScalarProduct(wannier_data, orbital_overlaps, 0, i, 0, j)<< " ";
    }
    std::cout << "\n";
  }

  std::cout << "\n\nWannier overlaps (different sites:)\n";

  for (size_t i  = 0; i < n_orbitals; ++i) {
    for (size_t j  = 0; j < n_orbitals; ++j) {
      std::cout << ScalarProduct(wannier_data, orbital_overlaps, 0, i, 1, j)<< " ";
    }
    std::cout << "\n";
  }
  


  std::cout << "\n\nTest finished...\n\n";

  return 0;
}
