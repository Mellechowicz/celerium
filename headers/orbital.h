#ifndef ORBITAL_H
#define ORBITAL_H

#include "arithmeticvector.h"

namespace celerium {

struct orbital {
  size_t n;
  size_t l;
  ArithmeticVector<3, double> position;
};

} // end namespace celerium

#endif /* ORBITAL_H */
