#ifndef BASIS_H
#define BASIS_H

#include "arithmeticvector.h"

namespace celerium {

class Basis {
 public:

  // Constructors.

  Basis() {}
  Basis(const std::array<ArithmeticVector, 3> &basis) {this->Set(basis);}

  // Operators.
  
  Basis &operator=(const Basis &rhs) {
    this->basis = rhs.basis;
    return *this;
  }

  virtual ~Basis() {};

  // Setters.
  
  void Set(const std::array<ArithmeticVector, 3> &basis) {
        double volume = fabs(basis[0]*(basis[1]^basis[2]));
    
    if (volume < 1e-10)
      throw std::invalid_argument("celerium::Basis: Basis \
vectors must be linearly independent.");
    this->basis = basis;
  }
  
  // Getters.

  const std::array<ArithmeticVector, 3> &GetVectors() const {return this->basis;}

 private:
  std::array<ArithmeticVector, 3> basis;
};

}  // end namspace celerium

#endif /* BASIS_H */
