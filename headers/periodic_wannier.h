#ifndef PERIODIC_WANNIER_H
#define PERIODIC_WANNIER_H

#include "gslmatrixcomplex.h"
#include "arithmeticvector.h"
#include <vector>
#include <algorithm>
#include <fstream>

namespace celerium {

// Structure describing the position of the unit cell
// relative to the system origin.
//  * basis_coords: coordinates of the cell
//    in the Bravais lattice basis, e.g., (0, 0, 0) corresponds
//    the unit cell int he origin.
//  * absolute position: position of the cell in absolute units
//    (angstroms).
struct cell_position_struct {
  std::array<int, 3> basis_coords;
  ArithmeticVector absolute_position;
};

// Structure holding the weight of an orbital in the Wannier function
// and additional information required to identify the orbital.
//   * coeff: value of the coefficient.
//   * position_index: index of the unit cell to which the orbital belongs.
//   * orbital_index: index of an orbital.
struct extended_coeff_struct {
  double coeff;
  size_t position_index;
  size_t orbital_index;
};

// Structure describing single Wannier function.
//   * wannier_index: index identifying the Wannier function.
//     By construction of the orthogonzalization algorithm
//     those indices are related to the original orbtial indices,
//     e.g. Wannier function with wannier_index = 3 will have the largest
//     overlap with orbtial of index 3.
//   * basis_coords & absolute_position: meaning the same
//     as in cell_position_struct.
//   * extended_coeffs: vector encoding the orbitals and weights contributing
//     to the Wannier function.
struct wannier_struct {
  size_t wannier_index;
  std::array<int, 3> basis_coords;
  ArithmeticVector absolute_position;
  std::vector<extended_coeff_struct> extended_coeffs;
};

// Low-level function for calculating Wannier functions. It is used internally
// by the class WannierData.
//   * orbital_overlaps: vector containing overlaps between initial
//     set of orbitals. The first argument of the pair is distance between
//     unit cells (in the lattice basis), whereas the second denotes the
//     square matrix of orbtial overlaps for this cell distance. For instance,
//     a distance [1, 0, 0] and a matrix should read
//     [<psi^0_[0, 0, 0]|psi^0_[1, 0, 0]>,  <psi^1_[0, 0, 0]|psi^0_[1, 0, 0]>]
//     [<psi^0_[0, 0, 0]|psi^1_[1, 0, 0]>,  <psi^1_[0, 0, 0]|psi^1_[1, 0, 0]>],
//     where superscript corresponds to orbital index and subscript to the
//     unit cell position.
//     IMPORTANT REMARK: note that the overlap matrices for opposite
//     cell distances (e.g., [1, 2, 0] and [-1, -2, 0]) are trivially
//     related by transposition. Only one representative of each such pair
//     should be present in orbital_overlaps. The function will internally
//     handle this symmetry. E.g., for a 2D system with non-zero
//     overlaps between nearest neighbors the matrices for [0, 0, 0],
//     [1, 0, 0], [0, 1, 0] should be included (neither [-1, 0, 0] nor
//     [0, -1, 0] is present).
//   * orbital_range: the maximum range of orbitals contributing to the
//     Wannier functions.
//   * wannier_positions: requested cell positions for which Wannier functions
//     will be calculated.
//   * cell_positions: this vector will be filled with positions of the unit 
//     cells contributing to the Wannier functions.
//   * wannier_coefficients: the Wannier functions will be written to this
//     entry. The format is the following:
//     wannier_coefficients[i][j][k][l] contains the coefficients for the
//     Wannier in i-th requested position (defined in wannier_positions) with
//     index j (i.e., it will be close to j-th original orbital). The index k
//     denotes the orbital position defined in wannier_positions, whereas l is
//     orbital index. We can thus write schematically:
//     W[i, j] = \sum_{k, l} wannier_coefficients[i][j][k][l] \psi^l_{cell_positions[k]}
int PeriodicOrthogonalization(
    std::vector<
       std::pair<std::array<int, 3>,gsl::Matrix>> &orbital_overlaps,
    const std::array<int, 3> &orbital_range,
    const std::vector<std::array<int, 3>> &wannier_positions,
    std::vector<std::array<int, 3>> &cell_positions,
    std::vector<std::vector<std::vector<std::vector<double>>>> &wannier_coefficients);

// This is a class handling orthogonalization intended for external use.
class WannierData {
 public:

  // Initializes WannierData. Most of the input data are the same as for the 
  // function PeriodicOrthogonalization. Additional parameters are:
  //  * basis: real-space basis of the lattice.
  //  * coefficient_cutoff: lower cutoff for the the Wannier coefficients.
  //    Orbitals with smaller weights will be disregarded. Setting a non-zero
  //    cutoff may substantially speed up the calculations.
  void Initialize(
      std::vector<std::pair<std::array<int, 3>, gsl::Matrix>> &orbital_overlaps,
      const std::array<int, 3> &orbital_range,
      const std::vector<std::array<int, 3>> &wannier_positions,
      const std::array<ArithmeticVector, 3> &basis,
      double coefficient_cutoff);

  // Saves Wannier data to file.
  int SaveToFile(const char *file_name);

  // Loads Wannier data to file.
  int LoadFromFile(const char *file_name);

  const std::vector<cell_position_struct> &GetOrbitalPositions() const {
    return this->cell_positions;
  }

  const std::vector<wannier_struct> &GetWanniers() const {
    return this->wanniers;
  }

  const std::array<int, 3> GetOrbitalRange() const {
    return this->orbital_range;
  };

  WannierData &operator=(const WannierData &rhs);
  
 private:
  std::vector<cell_position_struct> cell_positions;
  std::vector<wannier_struct> wanniers;
  std::array<int, 3> orbital_range;
  double coefficient_cutoff;
};

// Computes the scalar product between two Wannier functions.
double ScalarProduct(
  const WannierData &wannier_data,
  const std::vector<std::pair<std::array<int, 3>,
                              gsl::Matrix>> &orbital_overlaps,
  size_t position1_index, size_t wannier1_index,
  size_t position2_index, size_t wannier2_index);

} // end namespace celerium

#endif /* PERIODIC_WANNIER_H */
