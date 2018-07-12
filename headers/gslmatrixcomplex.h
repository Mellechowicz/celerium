#ifndef _GSLMATRIXCOMPLEX_H
#define _GSLMATRIXCOMPLEX_H

#include <iostream>
#include <memory>
#include <stdexcept>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <complex>
#include "gslmatrix.h"

namespace celerium{
namespace gsl{

class VectorComplex{
	gsl_vector_complex* myself;
	bool KILLME;

public:
	VectorComplex();
	VectorComplex(size_t n, bool INITIALIZE = false);
	VectorComplex(const VectorComplex& rhs);
	VectorComplex(gsl_vector_complex* v);
	VectorComplex(const gsl_vector_complex* v);
	virtual ~VectorComplex();

	void zero();
	int one(size_t i);

	friend std::ostream& operator<<(std::ostream& str, const VectorComplex& A);
	void set(gsl_vector_complex* v);

	void resize(size_t n, bool INITIALIZE = false);

	const gsl_vector_complex& operator()() const;
	gsl_vector_complex* operator()();
	int swap(VectorComplex& rhs);

static int swap(VectorComplex& lhs, VectorComplex& rhs){
	return gsl_vector_complex_swap(lhs(), rhs());
}

        std::complex<double> operator()(size_t i) const;
  //double& operator()(size_t i);

	size_t size() const;

	/*
	 * ARITHM. OPERATORS
	 */

	VectorComplex& operator+=(const VectorComplex& rhs);
	VectorComplex operator+(const VectorComplex& rhs) const;
	VectorComplex& operator-=(const VectorComplex& rhs);
	VectorComplex operator-(const VectorComplex& rhs) const;
        std::complex<double> operator*(const VectorComplex& rhs) const;
}; // end of class VectorComplex

class MatrixComplex{
	gsl_matrix_complex *myself;
	bool KILLME;

public:
	MatrixComplex();
	MatrixComplex(size_t n, bool INITIALIZE = false);
	MatrixComplex(size_t n, size_t m, bool INITIALIZE = false);
	MatrixComplex(const MatrixComplex& rhs);
  MatrixComplex(size_t n, const std::initializer_list<std::complex<double>> init);
	virtual ~MatrixComplex();
	void zero();
	void one();

	friend std::ostream& operator<<(std::ostream& str, const MatrixComplex& A);
	void resize(size_t n, bool INITIALIZE = false);

	const gsl_matrix_complex& operator()() const;
	gsl_matrix_complex* operator()();
  //double& operator()(size_t i, size_t j);
  double& real(size_t i, size_t j);
  double& imag(size_t i, size_t j);
  std::complex<double> operator()(size_t i, size_t j) const;
	int swap(MatrixComplex& rhs);
static int swap(MatrixComplex& lhs, MatrixComplex& rhs){
	return gsl_matrix_complex_swap(lhs(), rhs());
}

	size_t rowNumber() const;
	size_t columnNumber() const;
	VectorComplex row(size_t i);
	VectorComplex column(size_t i);
	void invert();
	int transpose();

	/*
	 * ARITHM. OPERATORS
	 */

	MatrixComplex& operator+=(const MatrixComplex& rhs);
	MatrixComplex operator+(const MatrixComplex& rhs) const;
	MatrixComplex& operator-=(const MatrixComplex& rhs);
	MatrixComplex operator-(const MatrixComplex& rhs) const;
	MatrixComplex operator*(const MatrixComplex& rhs) const;
	VectorComplex operator*(const VectorComplex& rhs) const;
	friend VectorComplex operator*(const VectorComplex& lhs, const MatrixComplex& rhs);
	MatrixComplex& operator*=(const MatrixComplex& rhs);
	MatrixComplex& operator=(const MatrixComplex& rhs);
	MatrixComplex& operator=(const VectorComplex& rhs);
        MatrixComplex& operator*(std::complex<double> scalar);

	/*
	 * EIGENPROBLEM
	 */

	void symmetricEigenProblem(MatrixComplex& eigenvectors, Vector& eigenvalues);
  
	/*
	 * Apply function
	 */

	template<typename T>
	MatrixComplex apply(T F){
		MatrixComplex beta;
		Vector lambda;
		symmetricEigenProblem(beta,lambda);
		MatrixComplex result;
	        result= beta;
		result.invert();

		for(size_t i=0; i<lambda.size(); ++i){
		  double l = F(lambda(i));
		  if(!std::isfinite(l))
		    throw std::domain_error("Function provided returned non-finite result: "+std::to_string(l)+" !");
	   	  for(size_t j=0; j<lambda.size(); ++j){
		    result.real(j,i) *= l;
                    result.imag(j,i) *= l;
                  }
                  
		}

		result *= beta;
		return result;
	}
	

}; // end of class MatrixComplex

} // end of namespace gsl
} // end of namespace celerium

#endif
