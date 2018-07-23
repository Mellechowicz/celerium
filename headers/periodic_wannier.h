#ifndef PERIODIC_WANNIER_H
#define PERIODIC_WANNIER_H

#include "gslmatrixcomplex.h"
#include "arithmeticvector.h"
#include <vector>
#include <algorithm>
#include <fstream>

namespace celerium {

struct orbital_position_struct {
  std::array<int, 3> cell_position;
  ArithmeticVector absolute_position;
};

struct extended_coeff_struct {
  double coeff;
  size_t position_index;
  size_t orbital_index;
};

struct wannier_struct {
  size_t wannier_index;
  std::array<int, 3> cell_position;
  ArithmeticVector absolute_position;
  std::vector<extended_coeff_struct> extended_coeffs;
};

int PeriodicOrthogonalization(
    std::vector<std::pair<std::array<int, 3>,
                          gsl::Matrix>> &orbital_overlaps,
    const std::array<int, 3> &wannier_range,
    const std::vector<std::array<int, 3>> &wannier_positions,
    std::vector<std::array<int, 3>> &orbital_positions,
    std::vector<std::vector<std::vector<std::vector<double>>>> &coefficients) {
  
  const size_t n_orbitals = orbital_overlaps.front().second.rowNumber();

  const size_t n_supercell_sites =
      (2*wannier_range[0]+1)*(2*wannier_range[1]+1)*(2*wannier_range[2]+1);

  const size_t n_positions = wannier_positions.size();
  
  const std::complex<double> I {0.0, 1.0};
  
  coefficients =
      std::vector<std::vector<std::vector<std::vector<double>>>>(
          n_positions,
          std::vector<std::vector<std::vector<double>>>(
              n_orbitals,
              std::vector<std::vector<double>>(
                  n_supercell_sites,
                  std::vector<double>(n_orbitals, 0.0))));
     
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
        
    for (auto orbital_overlap : orbital_overlaps) {
          
      double k0_times_r0 = orbital_overlap.first[0]*k[0];
      double k1_times_r1 = orbital_overlap.first[1]*k[1];
      double k2_times_r2 = orbital_overlap.first[2]*k[2];

      std::complex<double> coeff =
          std::exp( I*(k0_times_r0 + k1_times_r1 + k2_times_r2) );
       
      // Fourier transform of overlap matrix.

      k_space_overlap += gsl::MatrixComplex(orbital_overlap.second)*coeff;

      if (orbital_overlap.first[0] != 0 ||
          orbital_overlap.first[1] != 0 ||
          orbital_overlap.first[2] != 0) {
        auto orbital_overlap_transposed = orbital_overlap.second;
        orbital_overlap_transposed.transpose(); 
        k_space_overlap +=
            gsl::MatrixComplex(orbital_overlap_transposed)*conj(coeff);
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
                            [i]
                            [orbital_index_2] += tmp;

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

////////////////////////////////////////////////////////////////////




class WannierData {
 public:

  void Initialize(
      std::vector<std::pair<std::array<int, 3>, gsl::Matrix>> &orbital_overlaps,
      const std::array<int, 3> &orbital_range,
      const std::vector<std::array<int, 3>> &wannier_positions,
      const std::array<ArithmeticVector, 3> &basis,
      double coefficient_cutoff) {

    std::vector<std::array<int, 3>> bare_orbital_positions;
    std::vector<std::vector<std::vector<std::vector<double>>>> bare_coefficients;

    PeriodicOrthogonalization(orbital_overlaps,
                             orbital_range,
                             wannier_positions,
                             bare_orbital_positions,
                             bare_coefficients);    
    
    for (size_t i = 0; i < bare_coefficients.size(); ++i) {
      for (size_t j = 0; j < bare_coefficients[i].size(); ++j) {
        this->wanniers.push_back(wannier_struct());
        this->wanniers.back().wannier_index = j;
        this->wanniers.back().cell_position = wannier_positions[i];
        this->wanniers.back().absolute_position =
            basis[0]*(double)wannier_positions[i][0] +
            basis[1]*(double)wannier_positions[i][1] +
            basis[2]*(double)wannier_positions[i][2];
        for (size_t k = 0; k < bare_coefficients[i][j].size(); ++k) {
          for (size_t l = 0; l < bare_coefficients[i][j][k].size(); ++l) {
            if (fabs(bare_coefficients[i][j][k][l]) >= coefficient_cutoff) {

              auto located_orbital_position = 
                  std::find_if(this->orbital_positions.begin(),
                               this->orbital_positions.end(),
                               [&](const orbital_position_struct o) {
                                 auto c1 = o.cell_position[0] ==
                                           bare_orbital_positions[k][0];
                                 auto c2 = o.cell_position[1] ==
                                           bare_orbital_positions[k][1];
                                 auto c3 = o.cell_position[2] ==
                                           bare_orbital_positions[k][2];
                                 return c1 && c2 && c3;
                               });

              if (located_orbital_position == this->orbital_positions.end()) {
                orbital_position_struct tmp_orbital_position;
                tmp_orbital_position.cell_position = bare_orbital_positions[k];
                tmp_orbital_position.absolute_position =
                    basis[0]*(double)bare_orbital_positions[k][0] +
                    basis[1]*(double)bare_orbital_positions[k][1] +
                    basis[2]*(double)bare_orbital_positions[k][2];
                this->wanniers.back().extended_coeffs.push_back(
                    {bare_coefficients[i][j][k][l],this->orbital_positions.size()});
                this->orbital_positions.push_back(tmp_orbital_position);
              }
              else {
                this->wanniers.back().extended_coeffs.push_back(
                    { bare_coefficients[i][j][k][l],
                          (size_t)std::distance(this->orbital_positions.begin(),
                                                located_orbital_position),
                      l});
              }
            }
          }
        }
      }
    }

    
  }
      
  int SaveToFile(const char *file_name) {
    std::ofstream file(file_name, std::ofstream::trunc);
    
    file << orbital_positions.size() << " ";
    file << wanniers.size() << " ";

    for (size_t i = 0; i < orbital_positions.size(); ++i) {
      file << orbital_positions[i].cell_position[0] << " ";
      file << orbital_positions[i].cell_position[1] << " ";
      file << orbital_positions[i].cell_position[2] << " ";
      file << orbital_positions[i].absolute_position[0] << " ";
      file << orbital_positions[i].absolute_position[1] << " ";
      file << orbital_positions[i].absolute_position[2] << " ";
    }

    for (size_t i = 0; i < wanniers.size(); ++i) {
      file << wanniers[i].wannier_index << " ";
      file << wanniers[i].cell_position[0] << " ";
      file << wanniers[i].cell_position[1] << " ";
      file << wanniers[i].cell_position[2] << " ";
      file << wanniers[i].absolute_position[0] << " ";
      file << wanniers[i].absolute_position[1] << " ";
      file << wanniers[i].absolute_position[2] << " ";
      file << wanniers[i].extended_coeffs.size() << " ";
      for (size_t j  = 0; j < wanniers[i].extended_coeffs.size(); ++j) {
        file << wanniers[i].extended_coeffs[j].coeff << " ";
        file << wanniers[i].extended_coeffs[j].position_index << " ";
        file << wanniers[i].extended_coeffs[j].orbital_index << " ";
      }
    }
    
    file.close();
    return 0;
  }

  int LoadFromFile(const char *file_name) {
    std::ifstream file(file_name);

    size_t n_orbital_positions;
    size_t n_wanniers;
    
    file >> n_orbital_positions;
    file >> n_wanniers;

    this->orbital_positions.resize(n_orbital_positions);
    this->wanniers.resize(n_wanniers);

    for (size_t i = 0; i < n_orbital_positions; ++i) {
      file >> orbital_positions[i].cell_position[0];
      file >> orbital_positions[i].cell_position[1];
      file >> orbital_positions[i].cell_position[2];
      file >> orbital_positions[i].absolute_position[0];
      file >> orbital_positions[i].absolute_position[1];
      file >> orbital_positions[i].absolute_position[2];
    }

    for (size_t i = 0; i < n_wanniers; ++i) {
      size_t n_extended_coeffs;
      file >> wanniers[i].wannier_index;
      file >> wanniers[i].cell_position[0];
      file >> wanniers[i].cell_position[1];
      file >> wanniers[i].cell_position[2];
      file >> wanniers[i].absolute_position[0];
      file >> wanniers[i].absolute_position[1];
      file >> wanniers[i].absolute_position[2];
      file >> n_extended_coeffs;
      wanniers[i].extended_coeffs.resize(n_extended_coeffs);
      for (size_t j  = 0; j < n_extended_coeffs; ++j) {
        file >> wanniers[i].extended_coeffs[j].coeff;
        file >> wanniers[i].extended_coeffs[j].position_index;
        file >> wanniers[i].extended_coeffs[j].orbital_index;
      }
    }
    
    file.close();
    return 0;
  }

  const std::vector<orbital_position_struct> &GetOrbitalPositions() const
   {return this->orbital_positions;}
  const std::vector<wannier_struct> &GetWanniers() const
   {return this->wanniers;}

  WannierData &operator=(const WannierData &rhs) {
    this->orbital_positions = rhs.orbital_positions;
    this->wanniers = rhs.wanniers;
    return *this;
  }

  void Truncate(
      std::vector<std::pair<std::array<int, 3>,
                  std::vector<size_t>>> requested_wanniers);
  
 private:
  std::vector<orbital_position_struct> orbital_positions;
  std::vector<wannier_struct> wanniers;
  double cutoff;

  
};







} // end namespace celerium

#endif /* PERIODIC_WANNIER_H */
