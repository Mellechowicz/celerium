#include "../headers/lattice.h"




/*
void Lattice::AddElement(std::string name,
                         ArithmeticVector position_in_the_unit_cell,
                         std::vector<double> radial_potential_grid,
                         std::vector<double> radial_potential_values) {

  std::function name_search_function = [&](const internal_element_t &element)
                                       {return element.name == name;};
  
  auto located_element = std::find_if(this->elements.begin(),
                                      this->elements.end(),
                                      name_search_function);

  if (located_element != this->elements.end())
    throw std::invalid_argument("Lattice::AddElement: Attempted to \
add two elements with the same name.");

  internal_element_t new_internal_element;

  new_internal_element.name = name;
  new_internal_element.position_in_the_unit_cell = position_in_the_unit_cell;
  

  std::vector<sample_struct> samples;

  for (size_t i = 0; i < radial_potential_grid.size(); ++i) {
    samples.push_back({radial_potential_grid[i], radial_potential_values[i]});
  }

  std::function radial_potential_search_function =
      [&](const Interpolator &interp)
      {return interp.GetSamples() == samples;};
  
  auto located_radial_potential =
      std::find_if(this->radial_potentials.begin(),
                   this->radial_potentials.end(),
                   radial_potential_search_function);

  if (located_radial_potential == this->radial_potentials.end()) {
    Interpolator new_potential(samples);
    this->radial_potentials.emplace_back(new_potential);
    new_internal_element.radial_potential_index =
        this->radial_potentials.size() - 1;
  }
  else {
    new_internal_element.radial_potential_index =
        std::distance(this->radial_potentials.begin(),
                      located_radial_potential);    
  }

  this->elements.push_back(new_internal_element);
  

  
  
pa
}


*/


