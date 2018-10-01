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
        
    double volume =
        fabs( (basis_vectors[0]^basis_vectors[1])*
              basis_vectors[2] );
    
    double width0 =
        (basis_vectors[1]^basis_vectors[2]).length();
    
    double width1 =
        (basis_vectors[0]^basis_vectors[2]).length();
    
    double width2 =
        (basis_vectors[0]^basis_vectors[1]).length();
    
    width0 = fabs(volume/width0);
    width1 = fabs(volume/width1);
    width2 = fabs(volume/width2);
        
    int n0 = std::ceil(cutoff_radius/width0) + 1;
    int n1 = std::ceil(cutoff_radius/width1) + 1;
    int n2 = std::ceil(cutoff_radius/width2) + 1;
    
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
              
              dr = elementary_cell.GetOrbitalPositionInElementaryCell(w1) -
                         elementary_cell.GetOrbitalPositionInElementaryCell(w2) +
                         relative_cell_position;

              if (dr.length() > cutoff_radius) continue;

              dr = relative_cell_position;

              //std::cout << i0 << " " << i1 << " " << i2 << " " << w1 << " " << w2 << " ";
              //std::cout << elementary_cell.GetOrbitalPositionInElementaryCell(w1) << " ";
              //std::cout << elementary_cell.GetOrbitalPositionInElementaryCell(w2) << "\n";

              if (w1 == w2 && i0 == 0 && i1 == 0 && i2 ==0) {
                this->hamiltonian_terms.push_back({w1, w2, E0, relative_cell_distance, dr});
                this->hamiltonian_terms.push_back({w1, w2, U, relative_cell_distance, dr});
              }
              else {
                this->hamiltonian_terms.push_back({w1, w2, T, relative_cell_distance, dr});
                this->hamiltonian_terms.push_back({w1, w2, Up, relative_cell_distance, dr});
                this->hamiltonian_terms.push_back({w1, w2, J, relative_cell_distance, dr});
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
          b3 =
          {
            std::make_pair(0.0,integration_cutoff +
                           term.absolute_distance.length() / 2.0),
            std::make_pair(0.0,M_PI),
            std::make_pair(0.0,2.0*M_PI),
          };

      std::vector<double> resN (1), errN (1), pN (1);  
      steps = 0;

      std::function<int(const double *, double *)> integrand;
      integrand = [&](const double *xx, double *ff) {
                    double w1, w2, dw2, v;

                    auto position2 =
                    this->lattice_ptr->
                    GetElementaryCell().
                    GetOrbitalPosition(
                        this->hamiltonian_terms[term_index].wannier2_index);

                    const auto &dr =
                         this->hamiltonian_terms[term_index].absolute_distance;
                    
                    double xx2 [] = {
                      xx[0]*sin(xx[1])*cos(xx[2]) + position2[0] + 0*dr[0] / 2.0,
                      xx[0]*sin(xx[1])*sin(xx[2]) + position2[1] + 0*dr[1] / 2.0,
                      xx[0]*cos(xx[1]) + position2[2] + 0*dr[2] / 2.0
                    };

                    double xx1[] = {
                      xx2[0] - dr[0],
                      xx2[1] - dr[1],
                      xx2[2] - dr[2]
                    };

                    
                    w1 = this->lattice_ptr->EvaluateWannier(0, this->hamiltonian_terms[term_index].wannier1_index, xx1);
                    w2 = this->lattice_ptr->EvaluateWannier(0, this->hamiltonian_terms[term_index].wannier2_index, xx2);
                    dw2 = this->lattice_ptr->EvaluateLaplacian(0, this->hamiltonian_terms[term_index].wannier2_index, xx2);
                    v = this->lattice_ptr->EvaluateCrystalPotential(xx2);

                    ff[0] = (-3.80998208024*w1*dw2 + w1*w2*v)*xx[0]*xx[0]*sin(xx[1]);
                    return 0;
                  };

      engine.divonne_result(integrand, b3, resN, errN, pN, steps);

      value = resN[0];
      error = errN[0];
      return;
    }
    else if (term.type == U || term.type == Up) {
      
      std::vector<std::pair<double,double>>
          b6 =
      {
        std::make_pair(0.0,integration_cutoff),
        std::make_pair(0.0,M_PI),
        std::make_pair(0.0,2.0*M_PI),
        std::make_pair(0.0,integration_cutoff),
        std::make_pair(0.0,M_PI),
        std::make_pair(0.0,2.0*M_PI),
      };

      
      std::vector<double> resN (1), errN (1), pN (1);  
      steps = 0;

      std::function<int(const double *, double *)> integrand;
      integrand = [&](const double *xx, double *ff) {
                    double w1, w2;

                    auto position1 =
                    this->lattice_ptr->
                    GetElementaryCell().
                    GetOrbitalPosition(
                        this->hamiltonian_terms[term_index].wannier1_index);

                    auto position2 =
                    this->lattice_ptr->
                    GetElementaryCell().
                    GetOrbitalPosition(
                        this->hamiltonian_terms[term_index].wannier2_index);
                    
                    double xx1 [] = {
                      xx[0]*sin(xx[1])*cos(xx[2]) + position1[0],
                      xx[0]*sin(xx[1])*sin(xx[2]) + position1[1],
                      xx[0]*cos(xx[1]) + position1[2]
                    };

                    double xx2 [] = {
                      xx[3]*sin(xx[4])*cos(xx[5]) + position2[0],
                      xx[3]*sin(xx[4])*sin(xx[5]) + position2[1],
                      xx[3]*cos(xx[4]) + position2[2]
                    };

                    const auto &dr =
                         this->hamiltonian_terms[term_index].absolute_distance;
                    
                    double r = sqrt( (xx1[0]-xx2[0] + dr[0])*(xx1[0]-xx2[0] + dr[0]) +
                                     (xx1[1]-xx2[1] + dr[1])*(xx1[1]-xx2[1] + dr[1]) +
                                     (xx1[2]-xx2[2] + dr[2])*(xx1[2]-xx2[2] + dr[2]) );

                    if (r < 1e-10) r = 1e-10;                    
                    
                    w1 = this->lattice_ptr->EvaluateWannier(0, this->hamiltonian_terms[term_index].wannier1_index, xx1);
                    w2 = this->lattice_ptr->EvaluateWannier(0, this->hamiltonian_terms[term_index].wannier2_index, xx2);
                    ff[0] = 14.399645352 * w1*w1*w2*w2*xx[0]*xx[0]*sin(xx[1])*xx[3]*xx[3]*sin(xx[4]) / r;
                    return 0;
                  };

      engine.divonne_result(integrand, b6, resN, errN, pN, steps);

      value = resN[0];
      error = errN[0];
      return;      
    }
    else if (term.type == J) {
      
      std::vector<std::pair<double,double>>
          b6 =
      {
        std::make_pair(0.0,integration_cutoff),
        std::make_pair(0.0,M_PI),
        std::make_pair(0.0,2.0*M_PI),
        std::make_pair(0.0,integration_cutoff),
        std::make_pair(0.0,M_PI),
        std::make_pair(0.0,2.0*M_PI),
      };

      
      std::vector<double> resN (1), errN (1), pN (1);  
      steps = 0;

      std::function<int(const double *, double *)> integrand;
      integrand = [&](const double *xx, double *ff) {
        double w1x, w2x, w1xp, w2xp;

                    auto position1 =
                    this->lattice_ptr->
                    GetElementaryCell().
                    GetOrbitalPosition(
                        this->hamiltonian_terms[term_index].wannier1_index);

                    auto position2 =
                    this->lattice_ptr->
                    GetElementaryCell().
                    GetOrbitalPosition(
                        this->hamiltonian_terms[term_index].wannier2_index);
                    
                    double xx1 [] = {
                      xx[0]*sin(xx[1])*cos(xx[2]) + 0*position1[0],
                      xx[0]*sin(xx[1])*sin(xx[2]) + 0*position1[1],
                      xx[0]*cos(xx[1]) + 0*position1[2]
                    };

                    double xx2 [] = {
                      xx[3]*sin(xx[4])*cos(xx[5]) + 0*position2[0],
                      xx[3]*sin(xx[4])*sin(xx[5]) + 0*position2[1],
                      xx[3]*cos(xx[4]) + 0*position2[2]
                    };

                    const auto &dr =
                         this->hamiltonian_terms[term_index].absolute_distance;
                    
                    double r = sqrt( (xx1[0]-xx2[0] + dr[0])*(xx1[0]-xx2[0] + dr[0]) +
                                     (xx1[1]-xx2[1] + dr[1])*(xx1[1]-xx2[1] + dr[1]) +
                                     (xx1[2]-xx2[2] + dr[2])*(xx1[2]-xx2[2] + dr[2]) );

                    if (r < 1e-10) r = 1e-10;                    
                    
                    w1x = this->lattice_ptr->EvaluateWannier(0, this->hamiltonian_terms[term_index].wannier1_index, xx1);
                    w2x = this->lattice_ptr->EvaluateWannier(0, this->hamiltonian_terms[term_index].wannier2_index, xx1);
                    w1xp = this->lattice_ptr->EvaluateWannier(0, this->hamiltonian_terms[term_index].wannier1_index, xx2);
                    w2xp = this->lattice_ptr->EvaluateWannier(0, this->hamiltonian_terms[term_index].wannier2_index, xx2);

                    ff[0] = 2.0 * 14.399645352 * w1x * w1xp * w2x * w2xp *
                            xx[0]*xx[0]*sin(xx[1])*xx[3]*xx[3]*sin(xx[4]) / r;
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
