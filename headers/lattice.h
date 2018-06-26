#ifndef LATTICE_H
#define LATTICE_H

#include <orbital_class.h>
#include <element.h>

namespace celerium {


struct lattice_site_t {
  size_t element_index;
  ArithmeticVector position;
};


template <class Element>
class ElementaryCell {

 public:

  void SetBasis(const std::array<ArithmeticVector, 3> &basis) {

    double vol = fabs(basis[0]*(basis[1]^basis[2]));
    
    if (vol < 1e-7)
      throw std::invalid_argument("ElementaryCell::SetBasis: Basis \
vectors must be linearly independent.");

    this->basis = basis;
  }

  void AddSite(Element &element,
               ArithmeticVector position) {

    lattice_site_t lattice_site;
    lattice_site.position = position;
    
    auto located_element = 
        std::find_if(this->elements.begin(),
                     this->elements.end(),
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
    
    for (auto &orbital_class : element.GetOrbitalClasses()) {
      this->n_orbitals += orbital_class.GetActiveMValues().size();
    }

    this->lattice_sites.push_back(lattice_site);
    
  }

  size_t NOrbitals() const {return this->n_orbitals;}
 
  void EvaluateOrbitals(ArithmeticVector &coords,
                        std::vector<double> &result) {

    result.resize(this->n_orbitals);

    size_t i = 0;
    double radial_wf_value;

    for (auto &lattice_site : this->lattice_sites) {
      std::cout << "a\n";
      for (const auto &orbital_class :
               this->elements[lattice_site.element_index].GetOrbitalClasses()) {
        double r = (coords - lattice_site.position).length();
        radial_wf_value = orbital_class.GetRadialWF()(r);
        for (int m : orbital_class.GetActiveMValues()) {      
          result[i] = radial_wf_value *
                      RealSphericalHarmonic(orbital_class.GetL(),
                                            m,
                                            coords[0]-lattice_site.position[0],
                                            coords[1]-lattice_site.position[1],
                                            coords[2]-lattice_site.position[2]);
          i++;
        }
      }
    }
    
  }

  void EvaluateLaplacians(const ArithmeticVector &coords,
                          std::vector<double> result);
  

 private:
  size_t n_orbitals;
  std::vector<lattice_site_t> lattice_sites;
  std::vector<Element> elements;
  std::array<ArithmeticVector, 3> basis;
};











/*

struct site_t {
  std::string name;
  size_t element_index;
  ArithmeticVector position;
};

struct orbital_t {
  int l;
  int m;
  size_t site_index;
  size_t extended_radial_wf_index;
};


template <typename RadialPotential, typename RadialWaveFunction>
class Lattice {
 
 public:
  
  void SetRealSpaceBasis(std::vector<ArithmeticVector> real_space_basis) {


    
    if (real_space_basis.size() != 3)
      throw std::invalid_argument("Lattice::SetRealSpaceBasis: Invalid \
number of vector in the basis (must be 3).");


    double vol =
        fabs(real_space_basis[0]*(real_space_basis[1]^real_space_basis[2]));
    
    if (vol < 1e-7)
      throw std::invalid_argument("Lattice::SetRealSpaceBasis: Basis \
vectors must be linearly independent.");
    
    this->real_space_basis = real_space_basis;    
  }
  
  void AddElement(
      const Element<RadialPotential, RadialWaveFunction> &element) {

    if (element.GetName() == "")
      throw std::invalid_argument("Lattice::AddElement: Attempted to add \
element without a name.");
    
    auto name_search_function =
        [&](const Element<RadialPotential, RadialWaveFunction> &tested_element)
        {return tested_element.GetName() == element.GetName();};
  
    auto located_element = std::find_if(this->elements.begin(),
                                        this->elements.end(),
                                        name_search_function);

    if (located_element != this->elements.end())
      throw std::invalid_argument("Lattice::AddElement: Attempted to \
add two elements with the same name.");
    
    this->elements.push_back(element);
  }

  void AddSite(std::string site_name,
              std::string element_name,
              ArithmeticVector &site_position) {

    auto name_search_function =
        [&](const Element<RadialPotential, RadialWaveFunction> &tested_element)
        {return tested_element.GetName() == element_name;};
  
    auto located_element = std::find_if(this->elements.begin(),
                                        this->elements.end(),
                                        name_search_function);

    if (located_element == this->elements.end())
      throw std::invalid_argument("Lattice::AddSite: Attempted to \
add a lattice site with undefined element name.");

    auto site_name_search_function =
        [&](const site_t &tested_site)
        {return tested_site.name == site_name;};
  
    auto located_site = std::find_if(this->sites.begin(),
                                    this->sites.end(),
                                    site_name_search_function);

    if (located_site != this->sites.end())
      throw std::invalid_argument("Lattice::AddSite: Attempted to \
add two lattice sites with the same name.");

    site_t new_site;
    new_site.name = site_name;
    new_site.position = site_position;
    new_site.element_index = std::distance(this->elements.begin(),
                                                    located_element);

    this->sites.emplace_back(new_site);
  }

  void AddOrbital(std::string ion_name, int l, int m) {

    auto name_search_function =
        [&](const site_t &tested_element)
        {return tested_element.name == ion_name;};
  
    auto located_site = std::find_if(this->sites.begin(),
                                        this->sites.end(),
                                        name_search_function);

    if (located_site == this->sites.end())
      throw std::invalid_argument("Lattice::AddOrbital: Attempted to \
create orbital on latice site that does not exist (possible typo in name).");

    if (l < 0)
      throw std::invalid_argument("Lattice::AddOrbital: l must \
be non-negative.");

    if (abs(m) > l)
      throw std::invalid_argument("Lattice::AddOrbital: |m| \
must not exceed l");

    auto radial_search_function =
        [&](const extended_radial_wf<RadialWaveFunction> &tested_extended_wf)
        {return tested_extended_wf.l == l;};


    Element<RadialPotential, RadialWaveFunction> &
        ion_element = this->elements[located_site->element_index];
    
    auto located_wf = std::find_if(ion_element.GetRadialWaveFunctions().begin(),
                                   ion_element.GetRadialWaveFunctions().end(),
                                   radial_search_function);
    
    if (located_wf == ion_element.GetRadialWaveFunctions().end())
      throw std::invalid_argument("Lattice::AddOrbital: Attempted to \
add an orbital, for which the radial wave function was undefined.");

    orbital_t orbital;
    orbital.l = l;
    orbital.m = m;
    orbital.site_index =
        std::distance(this->sites.begin(),
                      located_site);
    orbital.extended_radial_wf_index =
        std::distance(ion_element.GetRadialWaveFunctions().begin(),
                      located_wf);
        
    orbitals.emplace_back(orbital);
  }

  size_t NOrbitals() const {
    return this->orbitals.size();
  }

  double WaveFunction(ArithmeticVector coords,
                      size_t orbital_index,
                      ArithmeticVectorN<3, int> cell) {

    auto &orbital = this->orbitals[orbital_index];
    const size_t ion_index = orbital.site_index;
    const size_t wf_index = orbital.extended_radial_wf_index;
    auto &ion = this->sites[ion_index];
    auto &element = this->elements[ion.element_index];
    auto &radial_wf = element.GetRadialWaveFunctions()[wf_index].radial_wf;

    const double r = coords.length();
    
    return radial_wf(r)*RealSphericalHarmonic(orbital.l, orbital.m,
                                 coords[0], coords[1],  coords[2]);
  }

  double CrystalPotential(ArithmeticVector coords, double cutoff_radius) {

    std::array<double, 3> proj;


    int n0 = abs(ceil(cutoff_radius / this->real_space_basis[0].length()));
    int n1 = abs(ceil(cutoff_radius / this->real_space_basis[1].length()));
    int n2 = abs(ceil(cutoff_radius / this->real_space_basis[2].length()));


    ArithmeticVector r_vec({0.0, 0.0, 0.0});
    
    double result = 0.0;
    for (const auto &orbital : orbitals) {
    
      const auto &ion =
          this->sites[orbital.site_index];

      for (size_t i = 0; i < 3; ++i) {
        proj[i] = ((coords-ion.position)*this->real_space_basis[i]) /
                  std::pow(this->real_space_basis[i].length(), 2);
        proj[i] = std::remainder(proj[i], 1.0);
      }
      
      for (int i0 = - n0; i0 <= n0; ++i0) {
        for (int i1 = - n1; i1 <= n1; ++i1) {
          for (int i2 = - n2; i2 <= n2; ++i2) {
 
            r_vec = (proj[0]+i0)*this->real_space_basis[0] +
                    (proj[1]+i1)*this->real_space_basis[1] +
                    (proj[2]+i2)*this->real_space_basis[2];

            //std::cout << r_vec.length() << "\n";
                
            if (r_vec.length() < cutoff_radius)            
              result += this->elements[ion.element_index].GetRadialPotential(r_vec.length());
          }
        }
      }
    }

    return result;

    
  }

  
  double Laplacian(ArithmeticVector coords,
                   size_t orbital_index,
                   ArithmeticVectorN<3, int> cell) {

    const auto &site =
        this->sites[this->orbitals[orbital_index].site_index];

    auto &element = this->elements[site.element_index];

    


    const size_t wf_index =
        this->orbitals[orbital_index].extended_radial_wf_index;

    ArithmeticVector r = coords - site.position -
                         this->real_space_basis[0] * cell[0] -
                         this->real_space_basis[1] * cell[1] -
                         this->real_space_basis[2] * cell[2];
                     
    
    double result = 0.0;
    result += element.GetRadialPotential(r.length());
    result -= element.GetRadialWaveFunctions()[wf_index].energy;
    result *= this->WaveFunction(r, orbital_index, cell);
    result *= 0.262468426082;  // h_bar^2 / 2 m_e
    
    return result;
    
  }
  
  
 private:
  std::vector<ArithmeticVector> real_space_basis;
  std::vector<Element<RadialPotential, RadialWaveFunction>> elements;
  std::vector<site_t> sites;
  std::vector<orbital_t> orbitals;
};


using LatticeI = Lattice<Interpolator, Interpolator>;


using LatticeF = Lattice<std::function<double(double)>,
                         std::function<double(double)>>;


*/

} // end namespace celerium



#endif /* LATTICE_H */
