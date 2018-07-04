#ifndef _GSLMATRIX_H
#define _GSLMATRIX_H

#include <iostream>
#include <memory>
#include <stdexcept>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

namespace celerium{
namespace gsl{

class Vector{
	gsl_vector* myself;
	bool KILLME;

public:
	Vector();
	Vector(size_t n, bool INITIALIZE = false);
	Vector(const Vector& rhs);
	Vector(gsl_vector* v);
	Vector(const gsl_vector* v);
	virtual ~Vector();

	void zero();
	int one(size_t i);

	friend std::ostream& operator<<(std::ostream& str, const Vector& A);
	void set(gsl_vector* v);

	void resize(size_t n, bool INITIALIZE = false);

	const gsl_vector& operator()() const;
	gsl_vector* operator()();
	int swap(Vector& rhs);

static int swap(Vector& lhs, Vector& rhs){
	return gsl_vector_swap(lhs(), rhs());
}

	double operator()(size_t i) const;
	double& operator()(size_t i);

	size_t size() const;

	/*
	 * ARITHM. OPERATORS
	 */

	Vector& operator+=(const Vector& rhs);
	Vector operator+(const Vector& rhs) const;
	Vector& operator-=(const Vector& rhs);
	Vector operator-(const Vector& rhs) const;
	double operator*(const Vector& rhs) const;
}; // end of class Vector

class Matrix{
	gsl_matrix *myself;
	bool KILLME;

public:
	Matrix();
	Matrix(size_t n, bool INITIALIZE = false);
	Matrix(size_t n, size_t m, bool INITIALIZE = false);
	Matrix(const Matrix& rhs);
	Matrix(size_t n, const std::initializer_list<double> init);
	virtual ~Matrix();
	void zero();
	void one();

	friend std::ostream& operator<<(std::ostream& str, const Matrix& A);
	void resize(size_t n, bool INITIALIZE = false);

	const gsl_matrix& operator()() const;
	gsl_matrix* operator()();
	double& operator()(size_t i, size_t j);
	double operator()(size_t i, size_t j) const;
	int swap(Matrix& rhs);
static int swap(Matrix& lhs, Matrix& rhs){
	return gsl_matrix_swap(lhs(), rhs());
}

	size_t rowNumber() const;
	size_t columnNumber() const;
	Vector row(size_t i);
	Vector column(size_t i);
	void invert();
	int transpose();

	/*
	 * ARITHM. OPERATORS
	 */

	Matrix& operator+=(const Matrix& rhs);
	Matrix operator+(const Matrix& rhs) const;
	Matrix& operator-=(const Matrix& rhs);
	Matrix operator-(const Matrix& rhs) const;
	Matrix operator*(const Matrix& rhs) const;
	Vector operator*(const Vector& rhs) const;
	friend Vector operator*(const Vector& lhs, const Matrix& rhs);
	Matrix& operator*=(const Matrix& rhs);
	Matrix& operator=(const Matrix& rhs);
	Matrix& operator=(const Vector& rhs);

	/*
	 * EIGENPROBLEM
	 */

	void symmetricEigenProblem(Matrix& eigenvectors, Vector& eigenvalues);

	/*
	 * Apply function
	 */

	template<typename T>
	Matrix apply(T F){
		Matrix beta;
		Vector lambda;
		symmetricEigenProblem(beta,lambda);
		Matrix result;
	        result= beta;
		result.invert();

		for(size_t i=0; i<lambda.size(); ++i){
		  double l = F(lambda(i));
		  if(!std::isfinite(l))
		    throw std::domain_error("Function provided returned non-finite result: "+std::to_string(l)+" !");
	   	  for(size_t j=0; j<lambda.size(); ++j)
		    result(j,i) *= l;
		}

		result *= beta;
		return result;
	}
	

}; // end of class Matrix

} // end of namespace gsl
} // end of namespace celerium

#endif
