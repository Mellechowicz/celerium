#ifndef PERIODIC_WANNIER_H
#define PERIODIC_WANNIER_H

#include "gslmatrixcomplex.h"
#include "arithmeticvector.h"
#include <vector>
#include <algorithm>

namespace celerium {

int PeriodicOthogonalization(
    std::vector<std::pair<std::array<int, 3>,
                          gsl::Matrix>> &real_space_overlaps,
    const std::array<int, 3> &wannier_range,
    const std::vector<std::array<int, 3>> &wannier_positions,
    std::vector<std::array<int, 3>> &orbital_positions,
    std::vector<std::vector<std::vector<double>>> &coefficients) {
  
  const size_t n_orbitals = real_space_overlaps.front().second.rowNumber();

  const size_t n_supercell_sites =
      (2*wannier_range[0]+1)*(2*wannier_range[1]+1)*(2*wannier_range[2]+1);

  const size_t n_positions = wannier_positions.size();
  
  const std::complex<double> I {0.0, 1.0};
  
  coefficients =
      std::vector<std::vector<std::vector<double>>>(
          n_positions,
          std::vector<std::vector<double>>(
              n_orbitals,
              std::vector<double>(n_supercell_sites*n_orbitals, 0.0)));
     
  gsl::MatrixComplex k_space_overlap(n_orbitals);
  gsl::MatrixComplex eigenvectors(n_orbitals);
  gsl::Vector eigenvalues(n_orbitals);
  
  std::vector<std::array<double, 3>> wave_vectors;

  for (int i0 = -(int)wannier_range[0]; i0 <= (int)wannier_range[0]; ++i0) {
    for (int i1 = -(int)wannier_range[1]; i1 <= (int)wannier_range[1]; ++i1) {
      for (int i2 = -(int)wannier_range[2]; i2 <= (int)wannier_range[2]; ++i2) {

        orbital_positions.push_back({i0, i1, i2});
        
        auto located_k =
            std::find_if(wave_vectors.begin(),
                         wave_vectors.end(),
                         [&](const std::array<double, 3> &k) {
                           auto c0 = (k[0] == -i0);
                           auto c1 = (k[1] == -i1);
                           auto c2 = (k[2] == -i2);
                           return c0 && c1 && c2;
                         });
        
        if (located_k != wave_vectors.end()) continue;
        
        wave_vectors.push_back({{(double)i0, (double)i1, (double)i2}});
      }
    }
  }

  for (auto &k : wave_vectors) {
    k[0] *= 2.0*M_PI/(2.0*wannier_range[0] + 1);
    k[1] *= 2.0*M_PI/(2.0*wannier_range[1] + 1);
    k[2] *= 2.0*M_PI/(2.0*wannier_range[2] + 1);
  }
  
  for (const auto &k : wave_vectors) {
      
    // Generating Fourier transforms.

    k_space_overlap.zero();
        
    for (auto real_space_overlap : real_space_overlaps) {
          
      double k0_times_r0 = real_space_overlap.first[0]*k[0];
      double k1_times_r1 = real_space_overlap.first[1]*k[1];
      double k2_times_r2 = real_space_overlap.first[2]*k[2];

      std::complex<double> coeff =
          std::exp( I*(k0_times_r0 + k1_times_r1 + k2_times_r2) );
       
      // Fourier transform of overlap matrix.

      k_space_overlap += gsl::MatrixComplex(real_space_overlap.second)*coeff;

      if (real_space_overlap.first[0] != 0 ||
          real_space_overlap.first[1] != 0 ||
          real_space_overlap.first[2] != 0) {
        auto real_space_overlap_transposed = real_space_overlap.second;
        real_space_overlap_transposed.transpose(); 
        k_space_overlap +=
            gsl::MatrixComplex(real_space_overlap_transposed)*conj(coeff);
      }
         
    }

    // We need to make a copy of overlaps as they will be
    // altered by eigensolver.
    auto k_space_overlap_copy = k_space_overlap;
    
    k_space_overlap_copy.symmetricEigenProblem(eigenvectors, eigenvalues);

    for (size_t i = 0; i < n_orbitals; ++i) {
      if (eigenvalues(i) <= 1e-14)
        throw std::runtime_error("celerium::PeriodicOrthogonalization: Provided overlap matrix is not positive definite.");
    }
        
    auto eigenvectors_transposed = eigenvectors;
    eigenvectors_transposed.hermitianConjugate(); 
       
    auto normalization_matrix = eigenvectors_transposed *
                                k_space_overlap *
                                eigenvectors;
        
    for (size_t i = 0; i < n_orbitals; ++i) {
      for (size_t j = 0; j < n_orbitals; ++j) {
        double norm = std::sqrt(abs(normalization_matrix(j, j)));
        eigenvectors.real(i, j) /= norm;
        eigenvectors.imag(i, j) /= norm;
      }
    }
        
    eigenvectors_transposed = eigenvectors;
    eigenvectors_transposed.hermitianConjugate();

    auto k_space_overlap_transposed = k_space_overlap;
    k_space_overlap_transposed.hermitianConjugate();

    size_t i {0};
    for (int i0 = -(int)wannier_range[0]; i0 <= (int)wannier_range[0]; ++i0) {
      for (int i1 = -(int)wannier_range[1]; i1 <= (int)wannier_range[1]; ++i1) {
        for (int i2 = -(int)wannier_range[2]; i2 <= (int)wannier_range[2]; ++i2) {
          
          auto entry =
              eigenvectors *
              eigenvectors_transposed *
              k_space_overlap_transposed *
              (k_space_overlap *
               eigenvectors *
               eigenvectors_transposed *
               k_space_overlap_transposed).apply([](double x) {
                   return 1.0/sqrt(x);
                 });

          // Naive wanniers
          // entry =
          //    eigenvectors;
          
          for (size_t pos_index = 0;
               pos_index < n_positions;
               ++pos_index) {

            double k0_times_r0 =
                (i0 - wannier_positions[pos_index][0]) * k[0];
            double k1_times_r1 =
                (i1 - wannier_positions[pos_index][1]) * k[1];
            double k2_times_r2 =
                (i2 - wannier_positions[pos_index][2]) * k[2];
            
            for (size_t orbital_index_1 = 0;
                 orbital_index_1 < n_orbitals;
                 ++orbital_index_1) {

              for (size_t orbital_index_2 = 0;
                   orbital_index_2 < n_orbitals;
                   ++orbital_index_2) {

                double tmp =
                    real(entry(orbital_index_2, orbital_index_1) *
                         std::exp(I*(k0_times_r0+k1_times_r1+k2_times_r2)) /
                         ((double)n_supercell_sites));

                if (k[0] != 0.0 || k[1] != 0.0 || k[2] != 0.0) tmp *= 2.0;

                coefficients[pos_index]
                            [orbital_index_1]
                            [orbital_index_2 + i*n_orbitals] += tmp;

              }
            }

          }
          ++i;
        }
      }
    }        
    
  }
 
  

  return 0;
}




////////////////////////////////////////////////////////////////////////////


double ScalarProduct(
    const std::vector<std::pair<std::array<int, 3>, gsl::Matrix>> &real_space_overlaps,
    const std::vector<std::array<int, 3>> &orbital_positions,
    const std::vector<std::vector<std::vector<double>>> &wanniers,
    size_t position1_index, size_t wannier1_index,
    size_t position2_index, size_t wannier2_index) {


  auto real_space_overlaps_extended {real_space_overlaps};
  for (const auto &real_space_overlap : real_space_overlaps) {
    if (real_space_overlap.first[0] != 0 ||
        real_space_overlap.first[1] != 0 ||
        real_space_overlap.first[2] != 0)
      real_space_overlaps_extended.push_back(
          { {-real_space_overlap.first[0],
             -real_space_overlap.first[1],
             -real_space_overlap.first[2]}, real_space_overlap.second});
    real_space_overlaps_extended.back().second.transpose();
  }
    
  int n0 {0};
  int n1 {0};
  int n2 {0};

  for (auto position : orbital_positions) {
    if (abs(position[0]) > n0) n0 = abs(position[0]);
    if (abs(position[1]) > n1) n1 = abs(position[1]);
    if (abs(position[2]) > n2) n2 = abs(position[2]);
  }

  int n_orbital_positions = (int)orbital_positions.size();
  
  size_t n_orbitals = real_space_overlaps[0].second.columnNumber();

  double result = 0;
  
  for (int  i = 0; i < n_orbital_positions; ++i) {
    for (const auto &real_space_overlap : real_space_overlaps_extended) {

      std::array<int, 3> r =
          { (orbital_positions[i][0] + real_space_overlap.first[0]),
            (orbital_positions[i][1] + real_space_overlap.first[1]),
            (orbital_positions[i][2] + real_space_overlap.first[2])};

        if (r[0] < -n0) r[0] += 2*n0 + 1;
        if (r[0] >  n0) r[0] -= 2*n0 + 1;
        if (r[1] < -n1) r[1] += 2*n1 + 1;
        if (r[1] >  n1) r[1] -= 2*n1 + 1;
        if (r[2] < -n2) r[2] += 2*n2 + 1;
        if (r[2] >  n2) r[2] -= 2*n2 + 1;

      int j = (r[0]+n0)*(2*n1+1)*(2*n2+1) + (r[1]+n1)*(2*n2+1) + (r[2]+n2);
      
      for (size_t o1 = 0; o1 < n_orbitals; ++o1) {
        for (size_t o2 = 0; o2 < n_orbitals; ++o2) {
          result +=
              wanniers[position1_index][wannier1_index][i*n_orbitals+o1]*
              wanniers[position2_index][wannier2_index][j*n_orbitals+o2]*
              real_space_overlap.second(o1, o2);       
          
        }
      }
      
    }
  }

  return result;
}

} // end namespace celerium

#endif /* PERIODIC_WANNIER_H */
