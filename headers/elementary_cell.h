#ifndef ELEMENTARY_CELL_H
#define ELEMENTARY_CELL_H

#include <basis.h>
#include <element.h>

namespace celerium {

// Structure internally representing
// the lattice site.
struct lattice_site_struct {
  std::string name;
  size_t element_index;
  ArithmeticVector position;
};


template <class Element = Element<Interpolator, OrbitalClass<Interpolator>>>
class ElementaryCell {

 public:

  ElementaryCell(const std::array<ArithmeticVector, 3> &basis) :
      n_orbitals(0),
      basis(basis) {}

  ElementaryCell(const ElementaryCell &e) :
      n_orbitals(e.n_orbitals),
      lattice_sites(e.lattice_sites),
      elements(e.elements),
      basis(e.basis) {}

  ElementaryCell &operator=(const ElementaryCell &rhs) {
    this->n_orbitals=rhs.n_orbitals;
    this->lattice_sites=rhs.lattice_sites;
    this->elements=rhs.elements;
    this->basis=rhs.basis;
  }
  
  void AddSite(const char *site_name,
               const Element &element,
               const ArithmeticVector &site_position) {

    auto located_site = std::find_if(this->lattice_sites.begin(),
                                     this->lattice_sites.end(),
                                     [&](lattice_site_struct ls) {
                                       return ls.name == site_name;
                                     });

    if (located_site != this->lattice_sites.end())
      throw std::invalid_argument("celerium::ElementaryCell::AddSite: \
Attempted to add two lattice sites with the same name.");
      
    lattice_site_struct lattice_site;
    lattice_site.name = site_name;
    lattice_site.position = site_position;
    
    auto located_element = 
        std::find_if(this->elements.begin(), this->elements.end(),
                     [&](Element &e) {
                       return e.GetName() == element.GetName();
                     });

    if (located_element == this->elements.end()) {
      this->elements.push_back(element);
      lattice_site.element_index = this->elements.size() - 1;
    }
    else {
      lattice_site.element_index =
          std::distance(this->elements.begin(), located_element);
    }
    
    for (auto &orbital_class : element.GetOrbitalClasses()) 
      this->n_orbitals += orbital_class.GetActiveMValues().size();
    
    this->lattice_sites.push_back(lattice_site);
  }

  size_t NOrbitals() const {return this->n_orbitals;}

  size_t NSites() const {return this->lattice_sites.size();}

  void EvaluateOrbitals(const ArithmeticVector &coords,
                        double result []) const {
 
    size_t i = 0;
    double radial_wf_value;

    for (const auto &lattice_site : this->lattice_sites) {
      double r = (coords - lattice_site.position).length();
      for (const auto &orbital_class :
               this->elements[lattice_site.element_index].GetOrbitalClasses()) {
        radial_wf_value = orbital_class.GetRadialWF()(r);
        for (int m : orbital_class.GetActiveMValues()) {
          result[i] = radial_wf_value *
                      RealSphericalHarmonic(orbital_class.GetL(),
                                            m,
                                            coords[0]-lattice_site.position[0],
                                            coords[1]-lattice_site.position[1],
                                            coords[2]-lattice_site.position[2]);
          i++;
        }  // active m values
      }  // orbital classes
    }  // lattice sites
  }

  void EvaluateOrbitals(const ArithmeticVector &coords,
                        std::vector<double> &result) const {
    result.resize(this->n_orbitals);
    EvaluateOrbitals(coords, result.data());
  }

  void EvaluateLaplacians(const ArithmeticVector &coords,
                          double result []) const {

    size_t i = 0;
    double radial_wf_value;
    double radial_potential_value;
    double energy;

    for (const auto &lattice_site : this->lattice_sites) {
      const Element &element = this->elements[lattice_site.element_index];
      double r = (coords - lattice_site.position).length();
      radial_potential_value = element.GetRadialPotential()(r);
      
      for (const auto &orbital_class :
               this->elements[lattice_site.element_index].GetOrbitalClasses()) {
        
        radial_wf_value = orbital_class.GetRadialWF()(r);
        energy = orbital_class.GetEnergy();
        
        for (int m : orbital_class.GetActiveMValues()) {
          result[i] = radial_potential_value;
          result[i] -= energy;
          result[i] *= radial_wf_value *
                       RealSphericalHarmonic(orbital_class.GetL(),
                                             m,
                                             coords[0]-lattice_site.position[0],
                                             coords[1]-lattice_site.position[1],
                                             coords[2]-lattice_site.position[2]);

          result[i] *= 0.262468426082;  // h_bar^2 / 2 m_e
          i++;
          
        } // active m values
      }  // orbital classes
    }  // lattice sites
  }

  void EvaluateLaplacians(const ArithmeticVector &coords,
                          std::vector<double> &result ) const {
      
    result.resize(this->n_orbitals);
    EvaluateLaplacians(coords, result.data());
  }


  // Not finished yet.
  /*
  double CellPotential(const ArithmeticVector &coords, double cutoff_radius) {
    
    std::array<double, 3> proj;

    int n0 = abs(ceil(cutoff_radius / this->basis.GetVectors()[0].length()));
    int n1 = abs(ceil(cutoff_radius / this->basis.GetVectors()[1].length()));
    int n2 = abs(ceil(cutoff_radius / this->basis.GetVectors()[2].length()));

    ArithmeticVector r_vec({0.0, 0.0, 0.0});
    
    double result = 0.0;
    for (const auto &lattice_site : this->lattice_sites) {

      for (size_t i = 0; i < 3; ++i) {
        proj[i] = ((coords-lattice_site.position)*this->basis[i]) /
                  std::pow(this->basis[i].length(), 2);
        proj[i] = std::remainder(proj[i], 1.0);
      }
      
      for (int i0 = - n0; i0 <= n0; ++i0) {
        for (int i1 = - n1; i1 <= n1; ++i1) {
          for (int i2 = - n2; i2 <= n2; ++i2) {
 
            r_vec = (proj[0]+i0)*this->basis[0] +
                    (proj[1]+i1)*this->basis[1] +
                    (proj[2]+i2)*this->basis[2];
            
            
            if (r_vec.length() < cutoff_radius) {
              //std::cout << "\nrvec: " << i0 << " " << i1 << " " << i2 << " " << r_vec.length()<< "\n";
              result += this->elements[lattice_site.element_index].GetRadialPotential()(r_vec.length());
              result -= this->elements[lattice_site.element_index].GetRadialPotential()(cutoff_radius);
            }
          }
        }
      }
    }

    return result;
  }
*/

 private:
  size_t n_orbitals; // Total number of orbitals in the unit cell.
  std::vector<lattice_site_struct> lattice_sites;  
  std::vector<Element> elements;
  Basis basis;
};

} // end namespace celerium

#endif /* ELEMENTARY_CELL_H */
