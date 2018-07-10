#ifndef LATTICE_H
#define LATTICE_H

template <class ElementaryCell>
Lattice {

  Lattice(ElementaryCell elementary_cell) :
      elementary_cell(elementary_cell) {}

  void Wannierize(std::array<int, 3> wannier_super_cell);


private:
  ElementaryCell elementary_cell;
  std::array<int, 3> wannier_super_cell
  
};


#endif /* LATTICE_H */
