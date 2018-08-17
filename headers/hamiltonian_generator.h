#ifndef HAMILTONIAN_GENERATOR_H
#define HAMILTONIAN_GENERATOR_H

#include "lattice.h"

namespace celerium {

enum hamiltonian_term_type {
  E0,
  T,
  U,
  Up,
  J,  
};

struct hamiltonian_term_struct {
  size_t wannier1_index;
  size_t wannier2_index;
  hamiltonian_term_type type;
  std::array<int, 3> cell_distance;
  ArithmeticVector absolute_distance;
};

template <class Lattice>
class HamiltonianGenerator {
 public:
  HamiltonianGenerator(double cutoff_radius, Lattice &lattice) {

    this->lattice_ptr = &lattice;

    this->cutoff_radius = cutoff_radius;

    const auto &elementary_cell = lattice.GetElementaryCell();
    const auto &basis_vectors = elementary_cell.GetBasis().GetVectors();
    
    int n0 =
        std::ceil(cutoff_radius/basis_vectors[0].length());
    int n1 =
        std::ceil(cutoff_radius/basis_vectors[1].length());
    int n2 =
        std::ceil(cutoff_radius/basis_vectors[2].length());
    
    size_t n_wanniers = lattice.NWanniers();
                            
    for (int i0 = -n0; i0 <= n0; ++i0) {
      for (int i1 = -n1; i1 <= n1; ++i1) {
        for (int i2 = -n2; i2 <= n2; ++i2) {
          auto relative_cell_position =
              basis_vectors[0]*(double)i0 +
              basis_vectors[1]*(double)i1 +
              basis_vectors[2]*(double)i2;

          std::array<int, 3> relative_cell_distance {i0, i1, i2};

          ArithmeticVector dr;
          for (size_t w1 = 0; w1 < n_wanniers; ++w1) {
            for (size_t w2 = 0; w2 < n_wanniers; ++w2) {
              
              dr = elementary_cell.GetSitePositionInElementaryCell(w1) -
                         elementary_cell.GetSitePositionInElementaryCell(w2) +
                         relative_cell_position;

              if (dr.length() > cutoff_radius) continue;

              if (w1 == w2 && i0 == 0 && i1 == 0 && i2 ==0) {
                this->hamiltonian_terms.push_back({w1, w2, E0, relative_cell_distance, dr});
                //this->hamiltonian_terms.push_back({w1, w2, U, relative_cell_distance, dr});
              }
              else {
                // this->hamiltonian_terms.push_back({w1, w2, T, relative_cell_distance, dr});
                //this->hamiltonian_terms.push_back({w1, w2, Up, relative_cell_distance, dr});
                //this->hamiltonian_terms.push_back({w1, w2, J, relative_cell_distance, dr});
              }
              
            }
          }
          
        }
      }
      
    }

    
  }

  void ComputeTerm(size_t term_index,
                   cuba::Cuba &engine,
                   double integration_cutoff,
                   double &value, double &error, int &steps) {

    const auto &term = this->hamiltonian_terms[term_index];
    
    if (term.type == T || term.type == E0) {

      std::vector<std::pair<double,double>>
          b3(3, std::make_pair(-integration_cutoff,integration_cutoff));
      std::vector<double> resN (1), errN (1), pN (1);  
      steps = 0;

      std::function<int(const double *, double *)> integrand;
      integrand = [&](const double *xx, double *ff) {
                    double w1, w2, dw2, v;
                    double xx_plus_dr []  = 
                        { xx[0] + this->hamiltonian_terms[term_index].absolute_distance[0],
                          xx[1] + this->hamiltonian_terms[term_index].absolute_distance[1],
                          xx[2] + this->hamiltonian_terms[term_index].absolute_distance[2] };

                    //ArithmeticVector xx_v({xx[0], xx[1], xx[2]});
                    //this->lattice_ptr->UpdateWanniers(xx_v);
                    //this->lattice_ptr->UpdateLaplacians(xx_v);

                    //w1 = this->lattice_ptr->GetWannier(0, this->hamiltonian_terms[term_index].wannier1_index);
                    //w2 = this->lattice_ptr->GetWannier(0, this->hamiltonian_terms[term_index].wannier2_index);
                    //dw2 = this->lattice_ptr->GetLaplacian(0, this->hamiltonian_terms[term_index].wannier2_index);
                    
                    w1 = this->lattice_ptr->EvaluateWannier(0, this->hamiltonian_terms[term_index].wannier1_index, xx_plus_dr);
                    w2 = this->lattice_ptr->EvaluateWannier(0, this->hamiltonian_terms[term_index].wannier2_index, xx);
                    dw2 = this->lattice_ptr->EvaluateLaplacian(0, this->hamiltonian_terms[term_index].wannier2_index, xx);
                    

                    v = this->lattice_ptr->EvaluateCrystalPotential(xx);
                    ff[0] = -3.80998208024*w1*dw2 + w1*w2*v;
                    return 0;
                  };

      engine.divonne_result(integrand, b3, resN, errN, pN, steps);

      value = resN[0];
      error = errN[0];
      return;
    }
    else if (term.type == U) {
      
      std::vector<std::pair<double,double>>
          b6(6, std::make_pair(-integration_cutoff,integration_cutoff));
      std::vector<double> resN (1), errN (1), pN (1);  
      steps = 0;

      std::function<int(const double *, double *)> integrand;
      integrand = [&](const double *xx, double *ff) {
                    double w1, w2, r;

                    r = std::sqrt( (xx[0] - xx[3])*(xx[0] - xx[3]) +
                                   (xx[1] - xx[4])*(xx[1] - xx[4]) +
                                   (xx[2] - xx[5])*(xx[2] - xx[5]) );
                    
                    if (r < 1e-10) r = 1e-10;
                    
                    w1 = this->lattice_ptr->EvaluateWannier(0, this->hamiltonian_terms[term_index].wannier1_index, xx);
                    w2 = this->lattice_ptr->EvaluateWannier(0, this->hamiltonian_terms[term_index].wannier2_index, xx + 3);
                    ff[0] = 0.5 * 0.5 * 14.399645352 * w1 * w1 * w2 * w2 / r;
                    return 0;
                  };

      engine.divonne_result(integrand, b6, resN, errN, pN, steps);

      value = resN[0];
      error = errN[0];
      return;
      
    }
    else if (term.type == J) {
 
      std::vector<std::pair<double,double>>
          b6(6, std::make_pair(-integration_cutoff,integration_cutoff));
      std::vector<double> resN (1), errN (1), pN (1);  
      steps = 0;

      std::function<int(const double *, double *)> integrand;
      integrand = [&](const double *xx, double *ff) {
                    double w1, w2, w3, w4, r;

                    r = std::sqrt( (xx[0] - xx[3])*(xx[0] - xx[3]) +
                                   (xx[1] - xx[4])*(xx[1] - xx[4]) +
                                   (xx[2] - xx[5])*(xx[2] - xx[5]) );

                    if (r < 1e-10) r = 1e-10;

                    double xx_plus_dr []  = 
                        { xx[0] + this->hamiltonian_terms[term_index].absolute_distance[0],
                          xx[1] + this->hamiltonian_terms[term_index].absolute_distance[1],
                          xx[2] + this->hamiltonian_terms[term_index].absolute_distance[2],
                          xx[3] + this->hamiltonian_terms[term_index].absolute_distance[0],
                          xx[4] + this->hamiltonian_terms[term_index].absolute_distance[1],
                          xx[5] + this->hamiltonian_terms[term_index].absolute_distance[2] };

                    w1 = this->lattice_ptr->EvaluateWannier(0, this->hamiltonian_terms[term_index].wannier1_index, xx_plus_dr);
                    w2 = this->lattice_ptr->EvaluateWannier(0, this->hamiltonian_terms[term_index].wannier2_index, xx + 3);
                    w3 = this->lattice_ptr->EvaluateWannier(0, this->hamiltonian_terms[term_index].wannier2_index, xx);
                    w4 = this->lattice_ptr->EvaluateWannier(0, this->hamiltonian_terms[term_index].wannier1_index, xx_plus_dr + 3);
                    ff[0] = 0.5 * 14.399645352 * w1 * w2 * w3 * w4 / r;
                    return 0;
                  };

      engine.divonne_result(integrand, b6, resN, errN, pN, steps);

      value = resN[0];
      error = errN[0];
      return;
    }
    else if (term.type == Up) {
 
      std::vector<std::pair<double,double>>
          b6(6, std::make_pair(-integration_cutoff,integration_cutoff));
      std::vector<double> resN (1), errN (1), pN (1);  
      steps = 0;

      std::function<int(const double *, double *)> integrand;
      integrand = [&](const double *xx, double *ff) {
                    double w1, w2, w3, w4, r;
                    
                    r = std::sqrt( (xx[0] - xx[3])*(xx[0] - xx[3]) +
                                   (xx[1] - xx[4])*(xx[1] - xx[4]) +
                                   (xx[2] - xx[5])*(xx[2] - xx[5]) );

                    if (r < 1e-10) r = 1e-10;

                    double xx_plus_dr []  = 
                        { xx[0] + this->hamiltonian_terms[term_index].absolute_distance[0],
                          xx[1] + this->hamiltonian_terms[term_index].absolute_distance[1],
                          xx[2] + this->hamiltonian_terms[term_index].absolute_distance[2],
                          xx[3] + this->hamiltonian_terms[term_index].absolute_distance[0],
                          xx[4] + this->hamiltonian_terms[term_index].absolute_distance[1],
                          xx[5] + this->hamiltonian_terms[term_index].absolute_distance[2] };
                    
                    if (r< 1e-10) r = 1e-10;
                    w1 = this->lattice_ptr->EvaluateWannier(0, this->hamiltonian_terms[term_index].wannier1_index, xx_plus_dr);
                    w2 = this->lattice_ptr->EvaluateWannier(0, this->hamiltonian_terms[term_index].wannier2_index, xx + 3);
                    w3 = this->lattice_ptr->EvaluateWannier(0, this->hamiltonian_terms[term_index].wannier1_index, xx_plus_dr);
                    w4 = this->lattice_ptr->EvaluateWannier(0, this->hamiltonian_terms[term_index].wannier2_index, xx + 3);
                    ff[0] = 0.5 * 0.5 * 14.399645352 * w1 * w2 * w3 * w4 / r;
                    return 0;
                  };

      engine.divonne_result(integrand, b6, resN, errN, pN, steps);

      value = resN[0];
      error = errN[0];
      return;
    }
    
  }
  
  size_t NTerms() {return this->hamiltonian_terms.size();}

  const std::vector<hamiltonian_term_struct> &GetTerms() {
    return this->hamiltonian_terms;
  }
  
 private:
  std::vector<hamiltonian_term_struct> hamiltonian_terms;
  double cutoff_radius;
  Lattice *lattice_ptr;
};

}  // celerium

#endif /* HAMILTONIAN_GENERATOR_H */
