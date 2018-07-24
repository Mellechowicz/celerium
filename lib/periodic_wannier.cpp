#include <periodic_wannier.h>

namespace celerium {

int PeriodicOrthogonalization(
    std::vector<std::pair<std::array<int, 3>,
    gsl::Matrix>> &orbital_overlaps,
    const std::array<int, 3> &wannier_range,
    const std::vector<std::array<int, 3>> &wannier_positions,
    std::vector<std::array<int, 3>> &cell_positions,
    std::vector<std::vector<std::vector<std::vector<double>>>> &wannier_coefficients) {
  
  const size_t n_orbitals = orbital_overlaps.front().second.rowNumber();

  const size_t n_supercell_sites =
      (2*wannier_range[0]+1)*(2*wannier_range[1]+1)*(2*wannier_range[2]+1);

  const size_t n_positions = wannier_positions.size();
  
  const std::complex<double> I {0.0, 1.0};
  
  wannier_coefficients =
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

        cell_positions.push_back({i0, i1, i2});
        
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
        throw std::runtime_error("PeriodicOrthogonalization: Provided overlap matrix is not positive definite.");
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

                wannier_coefficients[pos_index]
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

double ScalarProduct(
    const WannierData &wannier_data,
    const std::vector<std::pair<std::array<int, 3>,
    gsl::Matrix>> &orbital_overlaps,
    size_t position1_index, size_t wannier1_index,
    size_t position2_index, size_t wannier2_index) {

  double result {0};

  size_t n_orbitals = orbital_overlaps.front().second.columnNumber();
  
  const auto &wannier1 =
      wannier_data.GetWanniers()[position1_index*n_orbitals + wannier1_index];

  const auto &wannier2 =
      wannier_data.GetWanniers()[position2_index*n_orbitals + wannier2_index];

  for (const auto &wannier1_ext_coeff : wannier1.extended_coeffs) {
    for (const auto &wannier2_ext_coeff : wannier2.extended_coeffs) {
      const size_t &pos1_index =  wannier1_ext_coeff.position_index;
      const size_t &pos2_index =  wannier2_ext_coeff.position_index;
      const size_t &orb1_index =  wannier1_ext_coeff.orbital_index;
      const size_t &orb2_index =  wannier2_ext_coeff.orbital_index;
      const double &coeff1 =  wannier1_ext_coeff.coeff;
      const double &coeff2 =  wannier2_ext_coeff.coeff;

      const auto &pos1 =
          wannier_data.GetOrbitalPositions()[pos1_index].basis_coords;

      const auto &pos2 =
          wannier_data.GetOrbitalPositions()[pos2_index].basis_coords;

      std::array<int, 3> dr =
          {pos2[0] - pos1[0], pos2[1] - pos1[1], pos2[2] - pos1[2]};

      if (dr[0] > wannier_data.GetOrbitalRange()[0])
        dr[0] -= (2*wannier_data.GetOrbitalRange()[0] + 1);
      if (dr[0] < -wannier_data.GetOrbitalRange()[0])
        dr[0] += (2*wannier_data.GetOrbitalRange()[0] + 1);

      if (dr[1] > wannier_data.GetOrbitalRange()[1])
        dr[1] -= (2*wannier_data.GetOrbitalRange()[1] + 1);
      if (dr[1] < -wannier_data.GetOrbitalRange()[1])
        dr[1] += (2*wannier_data.GetOrbitalRange()[1] + 1);

      if (dr[2] > wannier_data.GetOrbitalRange()[2])
        dr[2] -= (2*wannier_data.GetOrbitalRange()[2] + 1);
      if (dr[2] < -wannier_data.GetOrbitalRange()[2])
        dr[2] += (2*wannier_data.GetOrbitalRange()[2] + 1);
      
      auto located_overlap = 
          std::find_if(orbital_overlaps.begin(),
                       orbital_overlaps.end(),
                       [&](const std::pair<std::array<int, 3>, gsl::Matrix> &o) {
                         return
                         o.first[0] == dr[0] &&
                         o.first[1] == dr[1] &&
                         o.first[2] == dr[2];
                       });

      auto located_overlap_inv = 
          std::find_if(orbital_overlaps.begin(),
                       orbital_overlaps.end(),
                       [&](const std::pair<std::array<int, 3>, gsl::Matrix> &o) {
                         return
                         o.first[0] == -dr[0] &&
                         o.first[1] == -dr[1] &&
                         o.first[2] == -dr[2] &&
                         ((dr[0] != 0) || (dr[1] != 0) || (dr[2] != 0));
                       });

      if (located_overlap != orbital_overlaps.end()) {
        result +=
            located_overlap->second(orb1_index, orb2_index)*coeff1*coeff2;
      }

      if (located_overlap_inv != orbital_overlaps.end()) {
        result +=
            located_overlap_inv->second(orb2_index, orb1_index)*coeff1*coeff2;
      }

    }
  }
  return result;
}


void WannierData::Initialize(
    std::vector<std::pair<std::array<int, 3>, gsl::Matrix>> &orbital_overlaps,
    const std::array<int, 3> &orbital_range,
    const std::vector<std::array<int, 3>> &wannier_positions,
    const std::array<ArithmeticVector, 3> &basis,
    double coefficient_cutoff) {

  std::vector<std::array<int, 3>> bare_cell_positions;
  std::vector<std::vector<std::vector<std::vector<double>>>> bare_coefficients;
  this->orbital_range = orbital_range;
  this->coefficient_cutoff = coefficient_cutoff;

  PeriodicOrthogonalization(orbital_overlaps,
                            orbital_range,
                            wannier_positions,
                            bare_cell_positions,
                            bare_coefficients);    
  
  for (size_t i = 0; i < bare_coefficients.size(); ++i) {
    for (size_t j = 0; j < bare_coefficients[i].size(); ++j) {
      this->wanniers.push_back(wannier_struct());
      this->wanniers.back().wannier_index = j;
      this->wanniers.back().basis_coords = wannier_positions[i];
      this->wanniers.back().absolute_position =
          basis[0]*(double)wannier_positions[i][0] +
          basis[1]*(double)wannier_positions[i][1] +
          basis[2]*(double)wannier_positions[i][2];
      for (size_t k = 0; k < bare_coefficients[i][j].size(); ++k) {
        for (size_t l = 0; l < bare_coefficients[i][j][k].size(); ++l) {
          if (fabs(bare_coefficients[i][j][k][l]) >= coefficient_cutoff) {

            auto located_cell_position = 
                std::find_if(this->cell_positions.begin(),
                             this->cell_positions.end(),
                             [&](const cell_position_struct o) {
                               auto c1 = o.basis_coords[0] ==
                               bare_cell_positions[k][0];
                               auto c2 = o.basis_coords[1] ==
                               bare_cell_positions[k][1];
                               auto c3 = o.basis_coords[2] ==
                               bare_cell_positions[k][2];
                               return c1 && c2 && c3;
                             });

            if (located_cell_position == this->cell_positions.end()) {
              cell_position_struct tmp_cell_position;
              tmp_cell_position.basis_coords = bare_cell_positions[k];
              tmp_cell_position.absolute_position =
                  basis[0]*(double)bare_cell_positions[k][0] +
                  basis[1]*(double)bare_cell_positions[k][1] +
                  basis[2]*(double)bare_cell_positions[k][2];
              this->wanniers.back().extended_coeffs.push_back(
                  {bare_coefficients[i][j][k][l],this->cell_positions.size()});
              this->cell_positions.push_back(tmp_cell_position);
            }
            else {
              this->wanniers.back().extended_coeffs.push_back(
                  { bare_coefficients[i][j][k][l],
                        (size_t)std::distance(this->cell_positions.begin(),
                                              located_cell_position),
                        l});
            }
          }
        }
      }
    }
  }

  
}

int WannierData::SaveToFile(const char *file_name) {
  std::ofstream file(file_name, std::ofstream::trunc);
  
  file << cell_positions.size() << " ";
  file << wanniers.size() << " ";
  file << orbital_range[0] << " ";
  file << orbital_range[1] << " ";
  file << orbital_range[2] << " ";
  file << coefficient_cutoff << " ";

  for (size_t i = 0; i < cell_positions.size(); ++i) {
    file << cell_positions[i].basis_coords[0] << " ";
    file << cell_positions[i].basis_coords[1] << " ";
    file << cell_positions[i].basis_coords[2] << " ";
    file << cell_positions[i].absolute_position[0] << " ";
    file << cell_positions[i].absolute_position[1] << " ";
    file << cell_positions[i].absolute_position[2] << " ";
  }

  for (size_t i = 0; i < wanniers.size(); ++i) {
    file << wanniers[i].wannier_index << " ";
    file << wanniers[i].basis_coords[0] << " ";
    file << wanniers[i].basis_coords[1] << " ";
    file << wanniers[i].basis_coords[2] << " ";
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

int WannierData::LoadFromFile(const char *file_name) {
  std::ifstream file(file_name);

  size_t n_cell_positions;
  size_t n_wanniers;
  
  file >> n_cell_positions;
  file >> n_wanniers;
  file >> orbital_range[0];
  file >> orbital_range[1];
  file >> orbital_range[2];
  file >> coefficient_cutoff;
  
  this->cell_positions.resize(n_cell_positions);
  this->wanniers.resize(n_wanniers);

  for (size_t i = 0; i < n_cell_positions; ++i) {
    file >> cell_positions[i].basis_coords[0];
    file >> cell_positions[i].basis_coords[1];
    file >> cell_positions[i].basis_coords[2];
    file >> cell_positions[i].absolute_position[0];
    file >> cell_positions[i].absolute_position[1];
    file >> cell_positions[i].absolute_position[2];
  }

  for (size_t i = 0; i < n_wanniers; ++i) {
    size_t n_extended_coeffs;
    file >> wanniers[i].wannier_index;
    file >> wanniers[i].basis_coords[0];
    file >> wanniers[i].basis_coords[1];
    file >> wanniers[i].basis_coords[2];
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

WannierData &WannierData::operator=(const WannierData &rhs) {
  this->cell_positions = rhs.cell_positions;
  this->wanniers = rhs.wanniers;
  this->orbital_range = rhs.orbital_range;
  this->coefficient_cutoff = rhs.coefficient_cutoff;
  return *this;
}




}  // celerium
