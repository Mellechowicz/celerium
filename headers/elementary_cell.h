#ifndef ELEMENTARY_CELL_H
#define ELEMENTARY_CELL_H

#include <basis.h>
#include <element.h>

namespace celerium {

// Structure internally representing
// the lattice site.
struct lattice_site_struct {
  std::string name;
  size_t element_index;         // Inequivalent elements are held in a separate container:
                                // this is the element index corresponding to the lattice site.
  ArithmeticVector position;    // Position of the site provided by the user.
  ArithmeticVector position_in_elementary_cell; // Guaranteed to be inside the
                                                // elementary cell. Used to
                                                // Efficiently calculate crystal
                                                // potential. Calculated automatically.
};

// Structure containing data used
// for computation of crystal potential.
struct crystal_potential_data_struct {
  std::vector<ArithmeticVector> contributing_cells; // Absolute Positions of
                                                    // cells contributing to the
                                                    // potential.
  double cutoff_radius; // Cutoff for pseudopotentials.
  double ionic_interacion_shift; // A constant shift of the potential
                                 // due to interaction between ions. 
};

// Structure containing data for quick access to orbital of given index.
struct orbital_handler_struct
{
  size_t lattice_site_index;
  size_t orbital_class_index;
  int m;
  std::string description;
};

template <class Element =
          Element<gsl::Interpolator, OrbitalClass<gsl::Interpolator>>>
class ElementaryCell {

 public:

  ElementaryCell(const std::array<ArithmeticVector, 3> &basis) :
      basis(basis) {}

  ElementaryCell(const ElementaryCell &e) :
      lattice_sites(e.lattice_sites),
      elements(e.elements),
      basis(e.basis),
      orbital_handlers(e.orbital_handlers),
      crystal_potential_data(e.crystal_potential_data) {}

  ElementaryCell &operator=(const ElementaryCell &rhs) {
    this->lattice_sites=rhs.lattice_sites;
    this->elements=rhs.elements;
    this->basis=rhs.basis;
    this->crystal_potential_data = rhs.crystal_potential_data;
    this->orbital_handlers = rhs.orbital_handlers;
  }

  void SetCrystalPotentialCutoff(double cutoff_radius) {

    this->crystal_potential_data.cutoff_radius = cutoff_radius;
   
    double volume =
        fabs( (this->basis.GetVectors()[0]^this->basis.GetVectors()[1])*
              this->basis.GetVectors()[2] );
    
    double width0 =
        (this->basis.GetVectors()[1]^this->basis.GetVectors()[2]).length();
    
    double width1 =
        (this->basis.GetVectors()[0]^this->basis.GetVectors()[2]).length();
    
    double width2 =
        (this->basis.GetVectors()[0]^this->basis.GetVectors()[1]).length();

    width0 = fabs(volume/width0);
    width1 = fabs(volume/width1);
    width2 = fabs(volume/width2);
        
    double n0 = std::ceil(cutoff_radius/width0) + 1;
    double n1 = std::ceil(cutoff_radius/width1) + 1;
    double n2 = std::ceil(cutoff_radius/width2) + 1;

    
    //std::cout << width0 << " " << width1 << " " << width2 << "\n";
    //std::cout << n0 << " " << n1 << " " << n2 << "\n";


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

    this->crystal_potential_data.ionic_interacion_shift = 0;
    
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
    /*
    for (size_t i = 0; i < 3; ++i) {
      double n =
          this->basis.GetVectors()[i]*lattice_site.position_in_elementary_cell;
      n /= this->basis.GetVectors()[i].length_squared();
      n = -std::floor(n);
      lattice_site.position_in_elementary_cell += this->basis.GetVectors()[i]*n;
    }
    */
    
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

    this->lattice_sites.push_back(lattice_site);
    
    size_t orbital_class_index = 0;
    orbital_handler_struct orbital_handler;
    for (const auto &orbital_class : element.GetOrbitalClasses()) {      
      orbital_handler.lattice_site_index = lattice_sites.size() - 1;
      orbital_handler.orbital_class_index = orbital_class_index;
      
      for (const auto &active_m_value : orbital_class.GetActiveMValues()) {
        orbital_handler.m = active_m_value;
        std::stringstream ss("");
        ss << lattice_site.name << "_"
           << orbital_class.GetN() <<
            OrbitalDescription(orbital_class.GetL(), active_m_value);
        orbital_handler.description = ss.str();
        this->orbital_handlers.push_back(orbital_handler);
      }
      orbital_class_index++;
    }
    
    
  }

  size_t NOrbitals() const {return this->orbital_handlers.size();}

  size_t NSites() const {return this->lattice_sites.size();}

  const std::string &GetOrbitalDescription(size_t orbital_index) {
    return this->orbital_handlers[orbital_index].description;
  }


  double EvaluateOrbital(size_t orbital_index,
                        const double coords []) const {
    const auto &orbital_handler =
        this->orbital_handlers[orbital_index];
    const auto &lattice_site =
        this->lattice_sites[orbital_handler.lattice_site_index];
    const auto &orbital_class =
        elements[lattice_site.element_index].
        GetOrbitalClasses()[orbital_handler.orbital_class_index];

    double x0 = coords[0] - lattice_site.position[0];
    double x1 = coords[1] - lattice_site.position[1];
    double x2 = coords[2] - lattice_site.position[2];
    double r = std::sqrt(x0*x0 + x1*x1 + x2*x2);

    return orbital_class.GetRadialWF()(r) *
        RealSphericalHarmonic(orbital_class.GetL(),
                              orbital_handler.m, x0, x1, x2);
  }

  double EvaluateLaplacian(size_t orbital_index,
                           const double coords []) const {
    const auto &orbital_handler = this->orbital_handlers[orbital_index];
    const auto &lattice_site =
        this->lattice_sites[orbital_handler.lattice_site_index];
    const auto &element = elements[lattice_site.element_index];
    const auto &orbital_class =
        element.GetOrbitalClasses()[orbital_handler.orbital_class_index];

    const double x0 = coords[0] - lattice_site.position[0];
    const double x1 = coords[1] - lattice_site.position[1];
    const double x2 = coords[2] - lattice_site.position[2];
    const double r = std::sqrt(x0*x0 + x1*x1 + x2*x2);

       
    double result = element.GetRadialPotential()(r);
    result -= orbital_class.GetEnergy();
    result *=
        orbital_class.GetRadialWF()(r) *
        RealSphericalHarmonic(orbital_class.GetL(),
                              orbital_handler.m, x0, x1, x2);

    result *= 0.262468426082;  // 1/(h_bar^2 / 2 m_e)

    return result;

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
    result.resize(this->NOrbitals());
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

  const ArithmeticVector &GetSitePosition(size_t site_index) const {
    return this->lattice_sites[site_index].position;
  }

  const ArithmeticVector &GetOrbitalPosition(size_t orbital_index) const {
    size_t site_index =  this->orbital_handlers[orbital_index].lattice_site_index;
    return lattice_sites[site_index].position;
  }


  const ArithmeticVector &GetSitePositionInElementaryCell(size_t site_index) const {
    return this->lattice_sites[site_index].position_in_elementary_cell;
  }

  const ArithmeticVector &GetOrbitalPositionInElementaryCell(size_t orbital_index) const {
    size_t site_index =  this->orbital_handlers[orbital_index].lattice_site_index;
    return lattice_sites[site_index].position_in_elementary_cell;
  }


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
  std::vector<lattice_site_struct> lattice_sites;  
  std::vector<Element> elements;
  Basis basis;
  std::vector<orbital_handler_struct> orbital_handlers;
  crystal_potential_data_struct crystal_potential_data;
};

} // end namespace celerium

#endif /* ELEMENTARY_CELL_H */

