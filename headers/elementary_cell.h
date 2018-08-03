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
  ArithmeticVector position;    // Position of the site provided by the user.
  ArithmeticVector position_in_elementary_cell; // Guaranteed to be inside the
                                                // elementary cell. Used to
                                                // efficiently calculate crystal
                                                // potential. 
};

// Structure containing data used
// for computation ofcrystal potential.
struct crystal_potential_data_struct {
  std::vector<ArithmeticVector> contributing_cells; // Absolute Positions of
                                                    // cells contributing to the
                                                    // potential.
  double cutoff_radius; // Cutoff for pseudopotentials.
  double ionic_interacion_shift;
};

template <class Element =
          Element<gsl::Interpolator, OrbitalClass<gsl::Interpolator>>>
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
    this->crystal_potential_data = rhs.crystal_potential_data;
  }

  void SetCrystalPotentialCutoff(double cutoff_radius) {

    this->crystal_potential_data.cutoff_radius = cutoff_radius;

    double n0 = std::ceil(cutoff_radius/this->basis.GetVectors()[0].length());
    double n1 = std::ceil(cutoff_radius/this->basis.GetVectors()[1].length());
    double n2 = std::ceil(cutoff_radius/this->basis.GetVectors()[2].length());

    for (double i0 = -n0; i0 <= n0; ++i0) {
      for (double i1 = -n1; i1 <= n1; ++i1) {
        for (double i2 = -n2; i2 <= n2; ++i2) {
          auto relative_position =
              this->basis.GetVectors()[0]*i0 +
              this->basis.GetVectors()[1]*i1 +
              this->basis.GetVectors()[2]*i2;
          
          this->crystal_potential_data.
              contributing_cells.push_back(relative_position);
        }
      }
    }

    this->crystal_potential_data.ionic_interacion_shift = 0.0;
    double r [] = {1.143561324124, 0.26543562, 0.312334532};
    this->crystal_potential_data.ionic_interacion_shift =
        -this->EvaluateCrystalPotential(r);
    
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
    lattice_site.position_in_elementary_cell = site_position;

    for (size_t i = 0; i < 3; ++i) {
      double n =
          this->basis.GetVectors()[i]*lattice_site.position_in_elementary_cell;
      n /= this->basis.GetVectors()[i].length_squared();
      n = -std::floor(n);
      lattice_site.position_in_elementary_cell += this->basis.GetVectors()[i]*n;
    }
    
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

  void GetOrbitalDescriptions(std::vector<std::string> &descriptions) {
    descriptions.clear();
    size_t i  = 0;
    for (const auto &lattice_site : this->lattice_sites) {
      for (const auto &orbital_class :
               this->elements[lattice_site.element_index].GetOrbitalClasses()) {
        for (int m : orbital_class.GetActiveMValues()) {
          std::stringstream ss("");
          ss << "index="<< i << ": " << lattice_site.name << ", "
             << orbital_class.GetN() <<
              OrbitalDescription(orbital_class.GetL(), m);
              descriptions.push_back(ss.str());
          ++i;
        }
      }
    }
  }


  
  void EvaluateOrbitals(const double coords [],
                        double result []) const {
 
    size_t i = 0;
    double radial_wf_value;
    double r, x0, x1, x2;
    
    for (const auto &lattice_site : this->lattice_sites) {

      x0 = coords[0] - lattice_site.position[0];
      x1 = coords[1] - lattice_site.position[1];
      x2 = coords[2] - lattice_site.position[2];  
      r = std::sqrt(x0*x0 + x1*x1 + x2*x2);
      
      for (const auto &orbital_class :
               this->elements[lattice_site.element_index].GetOrbitalClasses()) {
        radial_wf_value = orbital_class.GetRadialWF()(r);
        for (int m : orbital_class.GetActiveMValues()) {
          result[i] =
              radial_wf_value *
              RealSphericalHarmonic(orbital_class.GetL(), m, x0, x1, x2);
          i++;
        }  // active m values
      }  // orbital classes
    }  // lattice sites
  }

  void EvaluateOrbitals(const double coords [],
                        std::vector<double> &result) const {
    result.resize(this->n_orbitals);
    EvaluateOrbitals(coords, result.data());
  }

  void EvaluateLaplacians(const double coords [],
                          double result []) const {

    size_t i = 0;
    double radial_wf_value;
    double radial_potential_value;
    double energy;
    double r, x0, x1, x2;

    for (const auto &lattice_site : this->lattice_sites) {
      const Element &element = this->elements[lattice_site.element_index];

      x0 = coords[0] - lattice_site.position[0];
      x1 = coords[1] - lattice_site.position[1];
      x2 = coords[2] - lattice_site.position[2];  
      r = std::sqrt(x0*x0 + x1*x1 + x2*x2);
      
      radial_potential_value = element.GetRadialPotential()(r);
      
      for (const auto &orbital_class :
               this->elements[lattice_site.element_index].GetOrbitalClasses()) {
        
        radial_wf_value = orbital_class.GetRadialWF()(r);
        energy = orbital_class.GetEnergy();
        
        for (int m : orbital_class.GetActiveMValues()) {
          result[i] = radial_potential_value;
          result[i] -= energy;
          result[i] *=
              radial_wf_value *
              RealSphericalHarmonic(orbital_class.GetL(), m, x0, x1, x2);

          result[i] *= 0.262468426082;  // 1/(h_bar^2 / 2 m_e)
          i++;
          
        } // active m values
      }  // orbital classes
    }  // lattice sites
  }

  void EvaluateLaplacians(const double coords [],
                          std::vector<double> &result ) const {     
    result.resize(this->n_orbitals);
    EvaluateLaplacians(coords, result.data());
  }

  double EvaluateCrystalPotential(const double coords []) const {

    ArithmeticVector
        coords_in_elementary_cell({coords[0], coords[1], coords[2]});

    
    for (size_t i = 0; i < 3; ++i) {
      double n =
          this->basis.GetVectors()[i]*coords_in_elementary_cell;
      n /= this->basis.GetVectors()[i].length_squared();
      n = -std::floor(n);
      coords_in_elementary_cell += this->basis.GetVectors()[i]*n;
    }

    double result = this->crystal_potential_data.ionic_interacion_shift;
    ArithmeticVector relative_position;
    double dr;
    
    for (const auto &contributing_cell :
             this->crystal_potential_data.contributing_cells) {
      for (const auto &site : this->lattice_sites) {
        relative_position =
            coords_in_elementary_cell*(-1.0) +
            site.position_in_elementary_cell + contributing_cell;
        dr = relative_position.length();        
        if (dr < this->crystal_potential_data.cutoff_radius) 
          result += this->elements[site.element_index].GetRadialPotential()(dr);
      }
    }

    return result;
  }

  const Basis &GetBasis() const {return this->basis;}

 private:

  std::string OrbitalDescription(int l, int m) {
    std::string result;
    if (l == 0) {
      result = "s";
    }
    else if (l == 1) {
      switch (m) {
        case -1:
          result = "p_y";
          break;
        case 0:
          result = "p_z";
          break;
        case 1:
          result = "p_x";
          break;
      }
    }
    else if (l == 2) {
      switch (m) {
        case -2:
          result = "d_{xy}";
          break;
        case -1:
          result = "d_{yz}";
          break;
        case 0:
          result = "d_{z^2}";
          break;
        case 1:
          result = "d_{xz}";
          break;          
        case 2:
          result = "d_{x^2-y^2}";
          break;
      }
    }
    else if (l == 3) {
      switch (m) {
        case -3:
          result = "f_{y(3x^2-y^2)}";
          break;
        case -2:
          result = "f_{xyz}";
          break;
        case -1:
          result = "f_{yz^2}";
          break;
        case 0:
          result = "f_{z^3}";
          break;
        case 1:
          result = "f_{xz^2}";
          break;          
        case 2:
          result = "f_{z(x^2-y^2)}";
          break;
        case 3:
          result = "f_{x(x^2-3y^2)}";
          break;
      }
    }
    else {
      std::stringstream ss;
      ss << "l = " << l << ", m = " << m;
      result = ss.str();
    }
    return std::move(result);
}

  
  size_t n_orbitals; // Total number of orbitals in the unit cell.
  std::vector<lattice_site_struct> lattice_sites;  
  std::vector<Element> elements;
  Basis basis;
  crystal_potential_data_struct crystal_potential_data;
};

} // end namespace celerium

#endif /* ELEMENTARY_CELL_H */

