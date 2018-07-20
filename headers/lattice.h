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


struct wannier_data_struct {
  std::vector<std::array<int, 3>> wannier_positions;
  std::vector<std::array<int, 3>> orbital_positions;
  std::vector<ArithmeticVector> orbital_positions_absolute;
  std::vector<std::vector<std::vector<double>>> wannier_coefficients;

  int Save(const char* file_name) {
    std::ofstream file(file_name, std::ofstream::trunc);
    file << wannier_positions.size() << " ";
    file << orbital_positions.size() << " ";
    file << wannier_coefficients.front().size() << " \n\n";
    
    for (auto &position : wannier_positions)
      file << position[0] << " " << position[1] << " " << position[2] << " "; 
    file << "\n\n";
    for (auto &position : orbital_positions)
      file << position[0] << " " << position[1] << " " << position[2] << " ";
    file << "\n\n";
    for (auto &position : orbital_positions_absolute)
      file << position[0] << " " << position[1] << " " << position[2] << " ";
    file << "\n\n";
    for (const auto &coefficients_one_position : wannier_coefficients) {
      for (const auto &coefficients_one_wannier : coefficients_one_position) { 
        for (const auto &coefficient : coefficients_one_wannier) {
          file << coefficient << " ";
        }
        file << "\n\n";
      }
      file << "\n\n";
    }   
    file.close();
    return 0;
  }

  
  int Load(const char* file_name) {
    std::ifstream file(file_name);
    size_t n_wannier_pos;
    size_t n_orb_pos;
    size_t n_wannier;
    
    file >> n_wannier_pos >> n_orb_pos >> n_wannier;

    this->wannier_positions.resize(n_wannier_pos);
    this->orbital_positions.resize(n_orb_pos);
    this->orbital_positions_absolute.resize(n_orb_pos);
    this->wannier_coefficients =
        std::vector<std::vector<std::vector<double>>>(
            n_wannier_pos,
            std::vector<std::vector<double>>(
                n_wannier,
                std::vector<double>(n_orb_pos*n_wannier, 0.0)));

    for (size_t i_wannier_pos = 0; i_wannier_pos < n_wannier_pos; ++i_wannier_pos)
      file >> wannier_positions[i_wannier_pos][0] >>
              wannier_positions[i_wannier_pos][1] >>
              wannier_positions[i_wannier_pos][2]; 

    for (size_t i_orb_pos = 0; i_orb_pos < n_orb_pos; ++i_orb_pos)
      file >> orbital_positions[i_orb_pos][0] >>
              orbital_positions[i_orb_pos][1] >>
              orbital_positions[i_orb_pos][2];

    for (size_t i_orb_pos = 0; i_orb_pos < n_orb_pos; ++i_orb_pos)
      file >> orbital_positions_absolute[i_orb_pos][0] >>
              orbital_positions_absolute[i_orb_pos][1] >>
              orbital_positions_absolute[i_orb_pos][2]; 

    for (size_t i_wannier_pos = 0; i_wannier_pos < n_wannier_pos; ++i_wannier_pos) {
      for (size_t i_wannier = 0; i_wannier < n_wannier; ++i_wannier) { 
        for (size_t i = 0; i < n_wannier*n_orb_pos; ++i) {
          file >> wannier_coefficients[i_wannier_pos][i_wannier][i];
        }
        //file << "\n";
      }
      //file << "\n";
    }   
    
    file.close();
    return 0;
  }
  
};


template <class ElementaryCell>
class Lattice {

 public:

  Lattice(const ElementaryCell &elementary_cell) :
      elementary_cell(elementary_cell) {}

  int CalculateWannierCoefficients(
      const std::array<int, 3> &overlap_range,
      const std::array<int, 3> &wannier_range,
      const std::vector<std::array<int, 3>> &wannier_positions,
      cuba::Cuba &cuba_engine,
      bool verbose) {

    this->wannier_data.wannier_positions = wannier_positions;
    this->wannier_data.orbital_positions.clear();
    this->wannier_data.orbital_positions_absolute.clear();
    this->wannier_data.wannier_coefficients.clear();


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

          std::vector<std::pair<double,double>> b3(3,std::make_pair(-10,10));
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
      this->elementary_cell.EvaluateOrbitals(xx_arr, orbitals1);
      this->EvaluateOrbitals(xx_arr, drs[i_dr], orbitals2);
      
      size_t i = 0;
      for (size_t i1 = 0; i1 < n_orbitals; ++i1) {
        for (size_t i2 = 0; i2 < n_orbitals; ++i2) {
          ff[i] = orbitals1[i1]*orbitals2[i2];
          ++i;
        }
      }
      return 0;
    };


      cuba_engine.divonne_result(integrand, b3, resN, errN, pN, steps);

      if (verbose) {
        std::cerr << std::fixed << std::setprecision(4);
      
        std::cerr << "Position: (" << drs[i_dr][0] << ", " << drs[i_dr][1] << ", " <<
            drs[i_dr][2] << "). Overlap matrix " <<
            i_dr + 1 << " / " << drs.size() << ":\n";
      }
      
      
      
      for (size_t i  = 0; i < elementary_cell.NOrbitals(); ++i) {
        for (size_t j  = 0; j < elementary_cell.NOrbitals(); ++j) {
          overlap_matrix(i, j) = resN[i*elementary_cell.NOrbitals() + j];
          if (verbose)
            std::cerr << resN[i*elementary_cell.NOrbitals() + j] << " ";
        }
        if (verbose) std::cerr << "\n";
      }

      orbital_overlaps.push_back({drs[i_dr], overlap_matrix});

      if (verbose) std::cerr << "\n";
    }
    
    PeriodicOthogonalization(orbital_overlaps,
                             wannier_range,
                             this->wannier_data.wannier_positions,
                             this->wannier_data.orbital_positions,
                             this->wannier_data.wannier_coefficients);

    for (auto &orbital_position : this->wannier_data.orbital_positions) {
       auto r = 
        elementary_cell.GetBasis().GetVectors()[0]*(double)orbital_position[0]+
        elementary_cell.GetBasis().GetVectors()[1]*(double)orbital_position[1]+
        elementary_cell.GetBasis().GetVectors()[2]*(double)orbital_position[2];
      
      this->wannier_data.orbital_positions_absolute.push_back(r);
    }
    
    return 0;
  }

  void SaveWannierDataToFile(const char *file_name) {
    this->wannier_data.Save(file_name);
  }

  void LoadWannierDataFromFile(const char *file_name) {
    this->wannier_data.Load(file_name);
  }




  void EvaluateOrbitals(const ArithmeticVector &coords,
                        const std::array<int, 3> &orbital_position,
                        std::vector<double> &result) {

    auto &basis = elementary_cell.GetBasis().GetVectors();
    
    auto r = coords +
             basis[0]*((double)orbital_position[0])+
             basis[1]*((double)orbital_position[1])+
             basis[2]*((double)orbital_position[2]);
    
    this->elementary_cell.EvaluateOrbitals(r, result);
  }

  void UpdateWanniers(const ArithmeticVector &coords) {

    const size_t n_orbitals =
        this->elementary_cell.NOrbitals();
    const size_t n_orbital_positions =
        this->wannier_data.orbital_positions.size();

    this->evaluated_orbitals.resize(n_orbitals*n_orbital_positions);

    for (size_t i = 0; i < n_orbital_positions; ++i) {
      auto r = coords + this->wannier_data.orbital_positions_absolute[i];
      this->elementary_cell.EvaluateOrbitals(r,
                                             this->evaluated_orbitals.data()
                                             + i*n_orbitals);
    }
  }

  double GetWannier(size_t wannier_index, size_t wannier_position_index) {
    
    const auto &coeffs =
        this->wannier_data.wannier_coefficients[wannier_position_index]
                                               [wannier_index];

    double result = 0.0;
    for (size_t i = 0; i < coeffs.size(); ++i)
      result += coeffs[i]*this->evaluated_orbitals[i];

    return result;
  }
  
  void PrintWannier(size_t wannier_index, size_t wannier_position_index) {

    size_t n_orbitals = elementary_cell.NOrbitals();
  
    for (size_t orbital_position_index = 0;
         orbital_position_index < wannier_data.orbital_positions.size();
         ++orbital_position_index) {

      std::stringstream ss;
      ss << "(";
      ss << wannier_data.orbital_positions[orbital_position_index][0] << " ";
      ss << wannier_data.orbital_positions[orbital_position_index][1] << " ";
      ss << wannier_data.orbital_positions[orbital_position_index][2] << "):";
      
      std::cout << std::setw(8) << std::left << ss.str();
      
      for (size_t orbital_index = 0;
           orbital_index < n_orbitals;
           ++orbital_index) {
        std::cout << std::setw(6) <<
            wannier_data.wannier_coefficients[wannier_position_index][wannier_index]
            [orbital_index+n_orbitals*orbital_position_index] << " ";
      } // orbtial index
      std::cout << "\n";
    } // orbital positions
    std::cout << "\n";

  }

  
  private:
  ElementaryCell elementary_cell;
  wannier_data_struct wannier_data;
  std::vector<double> evaluated_orbitals;
};



}  // celerium
#endif /* LATTICE_H */
