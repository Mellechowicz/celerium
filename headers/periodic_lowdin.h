#ifndef PERIODIC_LOWDIN_H
#define PERIODIC_LOWDIN_H

#include "gslmatrix.h"
#include "arithmeticvector.h"
#include <vector>

/*
  std::array<ArithmeticVector, 3> k_space_basis;

  double cell_volume = basis[0] * (basis[1] ^ basis[2]);
  
  k_space_basis[0] = 2.0 * M_PI / cell_volume * basis[1] ^ basis[2];
  k_space_basis[1] = 2.0 * M_PI / cell_volume * basis[2] ^ basis[0];
  k_space_basis[2] = 2.0 * M_PI / cell_volume * basis[0] ^ basis[1];
*/


namespace celerium {

int PeriodicLowdin(std::vector<std::pair<std::array<int, 3>, gsl::Matrix>> &real_space_overlaps,
                   const std::array<size_t, 3> &supercell,
                   std::vector<std::vector<double>> &result) {

  size_t matrix_dim = real_space_overlaps.front().second.rowNumber();
  
  result =
      std::vector(real_space_overlaps.size(), std::vector<double>(matrix_dim, 0));

  gsl::Matrix k_space_overlap(matrix_dim);
  gsl::Matrix eigenvectors(matrix_dim);
  gsl::Vector eigenvalues(matrix_dim);

  for (size_t i_0 = 0; i_0 < supercell[0]; ++i_0) {
    for (size_t i_1 = 0; i_1 < supercell[1]; ++i_1) {
      for (size_t i_2 = 0; i_2 < supercell[2]; ++i_2) {

        k_space_overlap.zero();

        double k_0 = 2*M_PI*i_0/((double)supercell[0]);
        double k_1 = 2*M_PI*i_1/((double)supercell[1]);
        double k_2 = 2*M_PI*i_2/((double)supercell[2]);       
        
        for (auto real_space_overlap : real_space_overlaps) {
          
          const double k0_times_r0 = real_space_overlap.first[0]*k_0;
          const double k1_times_r1 = real_space_overlap.first[1]*k_1;
          const double k2_times_r2 = real_space_overlap.first[2]*k_2;

          double coeff = 1;
          coeff *= real_space_overlap.first[0] == 0 ? cos(k0_times_r0) : 2*cos(k0_times_r0);
          coeff *= real_space_overlap.first[1] == 0 ? cos(k1_times_r1) : 2*cos(k1_times_r1);
          coeff *= real_space_overlap.first[2] == 0 ? cos(k2_times_r2) : 2*cos(k2_times_r2);

          // Fourier transform of overlap matrix.
          k_space_overlap = real_space_overlap.second *
                            coeff + k_space_overlap;


        }

        k_space_overlap.symmetricEigenProblem(eigenvectors, eigenvalues);

        auto eigenvectors_transposed = eigenvectors;
        eigenvectors_transposed.transpose();
        auto k_space_overlap_transposed = k_space_overlap;
        k_space_overlap_transposed.transpose();
       
        auto normalization_matrix = eigenvectors_transposed *
                                    k_space_overlap_transposed *
                                    eigenvectors;


       
       for (size_t i = 0; i < matrix_dim; ++i) {
         for (size_t j = 0; j < matrix_dim; ++j) {
           double &element =  *gsl_matrix_ptr(eigenvectors(), i, j);
           element /= sqrt(normalization_matrix(j, j));
         }
       }
       std::cout << "norm: "<< (normalization_matrix(0, 0)) << "\n\n";

       std::cout << eigenvectors;
        
       eigenvectors_transposed = eigenvectors;
       eigenvectors_transposed.transpose();
      
       
        /*
        auto tmp_result = (k_space_overlap *
                       eigenvectors_transposed *
                       eigenvectors *
                       k_space_overlap_transposed);
                       
        
        tmp_result.apply([](double x) {return 1.0/std::sqrt(x);});
        */
        
        auto tmp_result = 
                     eigenvectors_transposed;


        //std::cout << k_space_overlap << "\n\n";

        size_t n_sites = supercell[0]*supercell[1]*supercell[2];
        
        for (size_t i = 0; i < real_space_overlaps.size(); ++i) {
          
          double k0_times_r0 =
              real_space_overlaps[i].first[0]*k_0;
          double k1_times_r1 =
              real_space_overlaps[i].first[1]*k_1;
          double k2_times_r2 =
              real_space_overlaps[i].first[2]*k_2;

          double entry = gsl_matrix_get(tmp_result(), 0, 0);



          std::cout << "entry: " << i  << " " << entry << " " << cos(k0_times_r0+k1_times_r1+k2_times_r2) << "\n\n";


          if (real_space_overlaps[i].first[0] == 0 &&
              real_space_overlaps[i].first[1] == 0 &&
              real_space_overlaps[i].first[2] == 0) {
          result[i][0] +=
              entry*cos(k0_times_r0+k1_times_r1+k2_times_r2)/n_sites;
          }
          else {
                      result[i][0] +=
              entry*cos(k0_times_r0+k1_times_r1+k2_times_r2)/n_sites*0.5;
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




} // end namespace celerium

#endif /* PERIODIC_LOWDIN_H */
