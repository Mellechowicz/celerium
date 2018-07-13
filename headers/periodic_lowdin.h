#ifndef PERIODIC_LOWDIN_H
#define PERIODIC_LOWDIN_H

#include "gslmatrixcomplex.h"
#include "arithmeticvector.h"
#include <vector>
#include <algorithm>

namespace celerium {



int PeriodicOthogonalization(
    std::vector<std::pair<std::array<int, 3>, gsl::MatrixComplex>> &real_space_overlaps,
    const std::array<size_t, 3> &supercell,
    const std::vector<std::array<int, 3>> &wannier_positions,
    std::vector<std::vector<std::vector<std::complex<double>>>> &result) {

  size_t matrix_dim = real_space_overlaps.front().second.rowNumber();

  const std::complex<double> I {0.0, 1.0};

  result =
      std::vector<std::vector<std::vector<std::complex<double>>>>(
          wannier_positions.size(),
          std::vector<std::vector<std::complex<double>>>(
              matrix_dim,
              std::vector<std::complex<double>>(
                  real_space_overlaps.size()*matrix_dim,
                  std::complex<double>(0.0))));

  gsl::MatrixComplex k_space_overlap(matrix_dim);
  gsl::MatrixComplex eigenvectors(matrix_dim);
  gsl::Vector eigenvalues(matrix_dim);
  std::array<double, 3> k;

  size_t i_0_max = supercell[0];
  size_t i_1_max = supercell[1];
  size_t i_2_max = supercell[2];

  /*  
  bool count_twice_k_0 ;
      
  if (i_0_max > 1) {
    i_0_max = i_0_max / 2 + i_0_max % 2;
    count_twice_k_0 = (i_0_max % 2 == 0 ? true : false);

    std::cout << "i0max: " << supercell[0] << "\n";
    std::cout << "i0max: " << i_0_max << "\n";
    
  }
  else if (i_1_max > 1) {
    i_1_max = i_1_max / 2 + i_1_max % 2;
    count_twice_k_0 = (i_1_max % 2 == 0 ? true : false);
  }
  else if (i_2_max > 1) {
    i_2_max = i_2_max / 2 + i_2_max % 2;
    count_twice_k_0 = (i_2_max % 2 == 0 ? true : false);
  }
  */
  
  for (size_t i_0 = 0; i_0 < i_0_max; ++i_0) {
    for (size_t i_1 = 0; i_1 < i_1_max; ++i_1) {
      for (size_t i_2 = 0; i_2 < i_2_max; ++i_2) {

        k[0] = 2*M_PI*i_0/((double)supercell[0]);
        k[1] = 2*M_PI*i_1/((double)supercell[1]);
        k[2] = 2*M_PI*i_2/((double)supercell[2]);               
      
        // Generating Fourier transforms.

        k_space_overlap.zero();
        
        for (auto real_space_overlap : real_space_overlaps) {
          
          double k0_times_r0 = real_space_overlap.first[0]*k[0];
          double k1_times_r1 = real_space_overlap.first[1]*k[1];
          double k2_times_r2 = real_space_overlap.first[2]*k[2];

          std::complex<double> coeff =
              std::exp( I*(k0_times_r0 + k1_times_r1 + k2_times_r2) );
       
          // Fourier transform of overlap matrix.
          k_space_overlap =
              k_space_overlap + real_space_overlap.second*coeff;
        }

        // We need to make a copy of overlaps as they will be
        // altered by eigensolver.
        auto k_space_overlap_copy = k_space_overlap;

        k_space_overlap_copy.symmetricEigenProblem(eigenvectors, eigenvalues);
        
        auto eigenvectors_transposed = eigenvectors;
        eigenvectors_transposed.hermitianConjugate(); 
       
        auto normalization_matrix = eigenvectors_transposed *
                                    k_space_overlap *
                                    eigenvectors;
        
       for (size_t i = 0; i < matrix_dim; ++i) {
         for (size_t j = 0; j < matrix_dim; ++j) {
           gsl_complex &element =
               *gsl_matrix_complex_ptr(eigenvectors(), i, j);
           double norm = std::sqrt(abs(normalization_matrix(j, j)));
           GSL_REAL(element) /= norm;
           GSL_IMAG(element) /= norm;
         }
       }
        
       eigenvectors_transposed = eigenvectors;
       eigenvectors_transposed.hermitianConjugate();

       auto k_space_overlap_transposed = k_space_overlap;
       k_space_overlap_transposed.hermitianConjugate();

       size_t n_sites = supercell[0]*supercell[1]*supercell[2];
        
        for (size_t i = 0; i < real_space_overlaps.size(); ++i) {
          
          auto entry =  eigenvectors *
                        eigenvectors_transposed *
                        k_space_overlap_transposed *
                        (k_space_overlap *
                         eigenvectors *
                         eigenvectors_transposed *
                         k_space_overlap_transposed).apply([](double x) {return 1.0/sqrt(x);});
          
          //auto entry = eigenvectors;



          for (size_t pos_index = 0;
               pos_index < wannier_positions.size();
               ++pos_index) {

            double k0_times_r0 =
                (real_space_overlaps[i].first[0] -
                 wannier_positions[pos_index][0]) * k[0];
            double k1_times_r1 =
                (real_space_overlaps[i].first[1] -
                 wannier_positions[pos_index][1]) * k[1];
            double k2_times_r2 =
                (real_space_overlaps[i].first[2] -
                 wannier_positions[pos_index][2]) * k[2];
            
            for (size_t orbital_index_1 = 0;
                 orbital_index_1 < matrix_dim;
                 ++orbital_index_1) {

              for (size_t orbital_index_2 = 0;
                   orbital_index_2 < matrix_dim;
                   ++orbital_index_2) {

                //if (k[0] == 0 && k[1] == 0 && k[2] == 0 && !count_twice_k_0) {
                
                   //result[pos_index][orbital_index_1][orbital_index_2 + i*matrix_dim] +=
                     // entry(orbital_index_2, orbital_index_1) *
                     //  std::cos(k0_times_r0+k1_times_r1+k2_times_r2) /
                     //   ((double)n_sites);
                 // }
                 //else {
                result[pos_index][orbital_index_1][orbital_index_2 + i*matrix_dim] +=
                    entry(orbital_index_2, orbital_index_1) *
                    std::exp(I*(k0_times_r0+k1_times_r1+k2_times_r2)) /
                    ((double)n_sites);                  
                //}
                
              }
              
            }
          
          }
          



      }        
    }
  }




    
}
  

  /*
  gsl::Matrix eigenvectors;
  gsl::Vector eigenvalues;

  for (size_t i = 0; i < overlap_matrices.size(); ++i) {
    overlap_matrices[i].second.symmetricEigenProblem(eigenvectors, eigenvalues);
  }
  */
  

  return 0;
}




////////////////////////////////////////////////////////////////////////////


  std::complex<double> overlap {0};


  size_t pos1 = 0;
  size_t pos2 = 0;
  size_t orb1 = 0;
  size_t orb2 = 0;


std::complex<double> ScalarProduct(
    const std::vector<std::pair<std::array<int, 3>, gsl::MatrixComplex>> &real_space_overlaps,
    const std::vector<std::vector<std::vector<std::complex<double>>>> &wanniers,
    size_t position1_index, size_t wannier1_index,
    size_t position2_index, size_t wannier2_index) {

  std::complex<double> overlap {0};
  
  int n = (int)real_space_overlaps.size();


  for (int  i = 0; i < n; ++i) {

    auto r_i = real_space_overlaps[i].first;
    
    for (int  j = 0; j < n; ++j) {
      
      auto r_j = real_space_overlaps[j].first;

      std::array<int, 3> dr = { r_j[0] - r_i[0],
                                r_j[1] - r_i[1],
                                r_j[2] - r_i[2] };

      auto located_overlap =
          std::find_if(real_space_overlaps.begin(),
                       real_space_overlaps.end(),
                       [&](const std::pair<std::array<int, 3>,
                           gsl::MatrixComplex> &arg) {
                         auto c0 = (arg.first[0] == dr[0] ||
                                    arg.first[0] == dr[0] + n ||
                                    arg.first[0] == dr[0] - n);
                         auto c1 = (arg.first[1] == dr[1] ||
                                    arg.first[1] == dr[1] + n||
                                    arg.first[1] == dr[1] - n);
                         auto c2 = (arg.first[2] == dr[2] ||
                                    arg.first[2] == dr[2] + n||
                                    arg.first[2] == dr[2] - n);
                         
                         return c0 && c1 && c2;
                       });

      if (located_overlap == real_space_overlaps.end()) continue;

      size_t n_orbitals = real_space_overlaps[0].second.columnNumber();
        
      for (size_t o1 = 0; o1 < n_orbitals; ++o1) {
        for (size_t o2 = 0; o2 < n_orbitals; ++o2) {     
          overlap +=
              conj(wanniers[position1_index][wannier1_index][i*n_orbitals+o1])*
              wanniers[position2_index][wannier2_index][j*n_orbitals+o2]*
              located_overlap->second(o1, o2);
        }
      }


      
    }
  }

  return overlap;
}

} // end namespace celerium

#endif /* PERIODIC_LOWDIN_H */
