#ifndef LATTICE_H
#define LATTICE_H

#include <vector>
#include <array>
#include <functional>
#include <periodic_wannier.h>
#include <newcubawrapper.h>
#include <elementary_cell.h>
#include <iomanip>
#include <fstream>

namespace celerium {

template <class ElementaryCell>
class Lattice {

 public:

  Lattice(const ElementaryCell &elementary_cell) :
      elementary_cell(elementary_cell) {}


  size_t NWanniers() const {
    return this->elementary_cell.NOrbitals();
  }

  const ElementaryCell &GetElementaryCell() const {
    return this->elementary_cell;
  }


  int CalculateWannierData(
      const std::array<int, 3> &overlap_range,
      const std::array<int, 3> &orbital_range,
      const std::vector<std::array<int, 3>> &wannier_positions,
      const std::vector<std::pair<double,double>> &integration_limits,
      cuba::Cuba &cuba_engine,
      double wannier_coefficient_cutoff,
      bool suppress_coefficiets_below_numerical_uncertainty,
      bool verbose) {

    std::vector<std::pair<std::array<int, 3>, gsl::Matrix>> orbital_overlaps;
    std::vector<std::array<int, 3>> drs;
    
    for (int i0 = -overlap_range[0]; i0 <= overlap_range[0]; ++i0) {
      for (int i1 = -overlap_range[1]; i1 <= overlap_range[1]; ++i1) {
        for (int i2 = -overlap_range[2]; i2 <= overlap_range[2]; ++i2) {   
          auto located_dr =
              std::find_if(drs.begin(), drs.end(),
                           [&](const std::array<int, 3> &dr) {
                             return i0 == -dr[0] &&
                                    i1 == -dr[1] &&
                                    i2 == -dr[2];});         
          if (located_dr != drs.end()) continue;
          drs.push_back({i0, i1, i2});          
        }
      }
    }
 
    for (size_t i_dr = 0; i_dr < drs.size(); ++i_dr) {
      size_t n_orbitals = elementary_cell.NOrbitals();

      std::vector<double> resN (n_orbitals*n_orbitals);
      std::vector<double> errN (n_orbitals*n_orbitals);
      std::vector<double> pN (n_orbitals*n_orbitals);  
      int steps = 0;
      gsl::Matrix overlap_matrix(n_orbitals);

      std::function<int(const double *, double *)> integrand;
      integrand = [&](const double *xx, double *ff) {
        std::vector<double> orbitals1 (n_orbitals);
        std::vector<double> orbitals2 (n_orbitals);
        ArithmeticVector xx_arr({xx[0], xx[1], xx[2]});
        this->elementary_cell.EvaluateOrbitals(xx, orbitals1);

        auto &basis = elementary_cell.GetBasis().GetVectors();
    
        auto r = xx_arr -
             basis[0]*((double)drs[i_dr][0])-
             basis[1]*((double)drs[i_dr][1])-
             basis[2]*((double)drs[i_dr][2]);

        this->elementary_cell.EvaluateOrbitals(r.data(), orbitals2);
      
        size_t i = 0;
        for (size_t i1 = 0; i1 < n_orbitals; ++i1) {
          for (size_t i2 = 0; i2 < n_orbitals; ++i2) {
            ff[i] = orbitals1[i1]*orbitals2[i2];
            ++i;
          }
        }
        return 0;
      };

      cuba_engine.divonne_result(integrand, integration_limits,
                                 resN, errN, pN, steps);

      if (verbose) {
        std::cerr << std::fixed << std::setprecision(4);
      
        std::cerr << "Position: (" << drs[i_dr][0]
                  << ", " << drs[i_dr][1] << ", " <<
            drs[i_dr][2] << "). Overlap matrix " <<
            i_dr + 1 << " / " << drs.size() << ":\n";
      }
      
      for (size_t i  = 0; i < elementary_cell.NOrbitals(); ++i) {
        for (size_t j  = 0; j < elementary_cell.NOrbitals(); ++j) {

          overlap_matrix(i, j) = resN[i*elementary_cell.NOrbitals() + j];

          if (suppress_coefficiets_below_numerical_uncertainty) {
            if (fabs(resN[i*elementary_cell.NOrbitals() + j]) <
                errN[i*elementary_cell.NOrbitals() + j])
              overlap_matrix(i, j) = 0;
          }
            
          if (verbose) {
            std::cerr << resN[i*elementary_cell.NOrbitals() + j] << " (";
            std::cerr << errN[i*elementary_cell.NOrbitals() + j] << ") ";
          }
        }
        if (verbose) std::cerr << "\n";
      }

      orbital_overlaps.push_back({drs[i_dr], overlap_matrix});

      if (verbose) std::cerr << "\n";
    }

    wannier_data.Initialize(orbital_overlaps,
                            orbital_range,
                            wannier_positions,
                            this->elementary_cell.GetBasis().GetVectors(),
                            wannier_coefficient_cutoff);

    size_t n_orbitals = this->elementary_cell.NOrbitals();
    size_t n_orbital_positions = this->wannier_data.GetOrbitalPositions().size();
    this->evaluated_orbitals.resize(n_orbitals*n_orbital_positions);
    this->evaluated_laplacians.resize(n_orbitals*n_orbital_positions);
        
    return 0;
  }

  
  void SaveWannierDataToFile(const char *file_name) {
    this->wannier_data.SaveToFile(file_name);
  }

  void LoadWannierDataFromFile(const char *file_name) {
    this->wannier_data.LoadFromFile(file_name);
    size_t n_orbtial_positions =
        this->wannier_data.GetOrbitalPositions().size();
    size_t n_orbitals = this->elementary_cell.NOrbitals();
    evaluated_orbitals.resize(n_orbitals*n_orbtial_positions);
    evaluated_laplacians.resize(n_orbitals*n_orbtial_positions);
  }

  void UpdateWanniers(const ArithmeticVector &coords) {

    const size_t n_orbitals = this->elementary_cell.NOrbitals();
    const size_t n_orbital_positions =
        this->wannier_data.GetOrbitalPositions().size();
      
    for (size_t i = 0; i < n_orbital_positions; ++i) {
      const auto r = coords -
                 this->wannier_data.GetOrbitalPositions()[i].absolute_position;
      this->elementary_cell.EvaluateOrbitals(r.data(),
                                             this->evaluated_orbitals.data() +
                                             i*n_orbitals);
    }
  }

  void UpdateLaplacians(const ArithmeticVector &coords) {

    const size_t n_orbitals = this->elementary_cell.NOrbitals();
    const size_t n_orbital_positions =
        this->wannier_data.GetOrbitalPositions().size();
      
    for (size_t i = 0; i < n_orbital_positions; ++i) {
      const auto r = coords -
                 this->wannier_data.GetOrbitalPositions()[i].absolute_position;
      this->elementary_cell.EvaluateLaplacians(r.data(),
                                               this->evaluated_laplacians.data() +
                                               i*n_orbitals);
    }
  }


  double GetWannier(size_t wannier_position_index, size_t wannier_index) const {

    const size_t n_orbitals = this->elementary_cell.NOrbitals();
    
    const auto &extended_coeffs =
        this->wannier_data.GetWanniers()
        [wannier_position_index*n_orbitals + wannier_index].extended_coeffs;

    double result {0};

    for (const auto& extended_coeff : extended_coeffs) {
      result += extended_coeff.coeff *
          this->evaluated_orbitals[n_orbitals*extended_coeff.position_index +
                                   extended_coeff.orbital_index];
    }
    
    return result;
  }

  void GetWanniers(size_t wannier_position_index, double result []) const {
    const size_t n_orbitals = this->elementary_cell.NOrbitals();
    const size_t shift = wannier_position_index*n_orbitals;
    
    for (size_t i = 0; i < n_orbitals; ++i) {
      const auto &extended_coeffs =
          this->wannier_data.GetWanniers()
          [shift + i].extended_coeffs;

      result[i] = 0;
      for (const auto& extended_coeff : extended_coeffs) {
        result[i] += extended_coeff.coeff *
             this->evaluated_orbitals[n_orbitals*extended_coeff.position_index +
                                      extended_coeff.orbital_index];
      }
    }
  }

  double EvaluateWannier(size_t wannier_position_index,
                       size_t wannier_index,
                       const double coords []) const {

    const size_t n_orbitals = this->elementary_cell.NOrbitals();
    
    const auto &extended_coeffs =
        this->wannier_data.GetWanniers()
        [wannier_position_index*n_orbitals + wannier_index].extended_coeffs;

   double result {0};


    for (const auto& extended_coeff : extended_coeffs) {   
      double xx [] = {
        -this->wannier_data.GetOrbitalPositions()[extended_coeff.position_index].absolute_position[0] + coords[0],
        -this->wannier_data.GetOrbitalPositions()[extended_coeff.position_index].absolute_position[1] + coords[1],
        -this->wannier_data.GetOrbitalPositions()[extended_coeff.position_index].absolute_position[2] + coords[2]
      };
      result += extended_coeff.coeff *
                this->elementary_cell.EvaluateOrbital(extended_coeff.orbital_index, xx);
    }
    
    return result;
  }

  
  double GetLaplacian(size_t wannier_position_index,
                      size_t wannier_index) const {

    const size_t n_orbitals = this->elementary_cell.NOrbitals();
    
    const auto &extended_coeffs =
        this->wannier_data.GetWanniers()
        [wannier_position_index*n_orbitals + wannier_index].extended_coeffs;

    double result {0};

    for (const auto& extended_coeff : extended_coeffs) {
      result += extended_coeff.coeff *
                this->evaluated_laplacians[
                    n_orbitals*extended_coeff.position_index +
                    extended_coeff.orbital_index];
    }
    
    return result;
  }

  void GetLaplacians(size_t wannier_position_index,
                     double result []) const {
    
    const size_t n_orbitals = this->elementary_cell.NOrbitals();
    const size_t shift = wannier_position_index*n_orbitals;
    
    for (size_t i = 0; i < n_orbitals; ++i) {
      const auto &extended_coeffs =
          this->wannier_data.GetWanniers()
          [shift + i].extended_coeffs;

      result[i] = 0;
      for (const auto& extended_coeff : extended_coeffs) {
        result[i] += extended_coeff.coeff *
                     this->evaluated_laplacians[
                         n_orbitals*extended_coeff.position_index +
                         extended_coeff.orbital_index];
      }
    }
  }  

  double EvaluateLaplacian(size_t wannier_position_index,
                       size_t wannier_index,
                       const double coords []) const {

    const size_t n_orbitals = this->elementary_cell.NOrbitals();
    
    const auto &extended_coeffs =
        this->wannier_data.GetWanniers()
        [wannier_position_index*n_orbitals + wannier_index].extended_coeffs;

   double result {0};

    for (const auto& extended_coeff : extended_coeffs) {
      double xx [] = {
        -this->wannier_data.GetOrbitalPositions()[extended_coeff.position_index].absolute_position[0] + coords[0],
        -this->wannier_data.GetOrbitalPositions()[extended_coeff.position_index].absolute_position[1] + coords[1],
        -this->wannier_data.GetOrbitalPositions()[extended_coeff.position_index].absolute_position[2] + coords[2]
      };
      result += extended_coeff.coeff *
                this->elementary_cell.EvaluateLaplacian(extended_coeff.orbital_index, xx);
    }
    
    return result;
  }



  double EvaluateCrystalPotential(const double coords []) const {
    return this->elementary_cell.EvaluateCrystalPotential(coords);
  }
  
  private:
  ElementaryCell elementary_cell;
  WannierData wannier_data;
  std::vector<double> evaluated_orbitals;
  std::vector<double> evaluated_laplacians;
};



}  // celerium

#endif /* LATTICE_H */
