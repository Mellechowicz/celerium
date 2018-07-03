#include "../headers/gslmatrix.h"

celerium::gsl::Vector::Vector(){
	myself = nullptr;
	KILLME = false;
}

celerium::gsl::Vector::Vector(size_t n, bool INITIALIZE){
	KILLME = true;
	if(INITIALIZE)  myself = gsl_vector_calloc(n);
	else		myself = gsl_vector_alloc(n);
}

celerium::gsl::Vector::Vector(const Vector& rhs){
	myself = gsl_vector_alloc(rhs().size);
	gsl_vector_memcpy(myself,&rhs());
	KILLME = true;
}

celerium::gsl::Vector::Vector(gsl_vector* v){
	KILLME = true;
	myself = v;
}

namespace celerium{ namespace gsl{
std::ostream& operator<<(std::ostream& str, const celerium::gsl::Vector& A){
	str<<"[ ";
	for(size_t i=0; i<A.size(); ++i)
	  str<<A(i)<<" ";
	str<<"]\n";

	return str;
}
}}

void celerium::gsl::Vector::set(gsl_vector* v){
	if(KILLME) gsl_vector_free(myself);
	myself = v;
	KILLME = true;
}

celerium::gsl::Vector::~Vector(){
	if(KILLME) gsl_vector_free(myself);
}

void celerium::gsl::Vector::resize(size_t n, bool INITIALIZE){
	if(KILLME)	gsl_vector_free(myself);
	if(INITIALIZE)  myself = gsl_vector_calloc(n);
	else		myself = gsl_vector_alloc(n);
	KILLME = true;
}

const gsl_vector& celerium::gsl::Vector::operator()() const{
	return *myself;
}

gsl_vector* celerium::gsl::Vector::operator()(){
	return myself; //unsafe!
}

int celerium::gsl::Vector::swap(celerium::gsl::Vector& rhs){
	return gsl_vector_swap(myself, rhs());
}

double celerium::gsl::Vector::operator()(size_t i) const{
	return gsl_vector_get(myself,i);
}

double& celerium::gsl::Vector::operator()(size_t i){
	return *gsl_vector_ptr(myself,i);
}

void celerium::gsl::Vector::zero(){
	gsl_vector_set_zero(myself);
}

int celerium::gsl::Vector::one(size_t i){
	return gsl_vector_set_basis(myself,i);
}

size_t celerium::gsl::Vector::size() const{
	return myself->size;
}

/*
 * ARITHM. OPERATORS
 */

celerium::gsl::Vector& celerium::gsl::Vector::operator+=(const celerium::gsl::Vector& rhs){
	if(size() != rhs.size())
	  throw std::logic_error("Oy! You can't add different-type vectors!");
	gsl_vector_add(myself,&rhs());
	return *this;
}

celerium::gsl::Vector celerium::gsl::Vector::operator+(const celerium::gsl::Vector& rhs) const{
	Vector result(*this);
	return result+=rhs;
}

celerium::gsl::Vector& celerium::gsl::Vector::operator-=(const celerium::gsl::Vector& rhs){
	if(size() != rhs.size())
	  throw std::logic_error("Oy! You can't substract different-type vectors!");
	gsl_vector_sub(myself,&rhs());
	return *this;
}

celerium::gsl::Vector celerium::gsl::Vector::operator-(const celerium::gsl::Vector& rhs) const{
	Vector result(*this);
	return result-=rhs;
}

double celerium::gsl::Vector::operator*(const celerium::gsl::Vector& rhs) const{
	double result;
	gsl_blas_ddot(myself,&rhs(),&result);
	return result;
}


celerium::gsl::Matrix::Matrix(){
	myself = nullptr;
	KILLME = false;
}

celerium::gsl::Matrix::Matrix(size_t n, bool INITIALIZE){
	  if(INITIALIZE) myself = gsl_matrix_calloc(n,n);
	  else		 myself = gsl_matrix_alloc(n,n);
	  KILLME = true;
}

celerium::gsl::Matrix::Matrix(size_t n, size_t m, bool INITIALIZE){
	if(INITIALIZE)  myself = gsl_matrix_calloc(n,m);
	else		myself = gsl_matrix_alloc(n,m);
	KILLME = true;
}

celerium::gsl::Matrix::Matrix(const celerium::gsl::Matrix& rhs){
	myself = gsl_matrix_alloc(rhs().size1,rhs().size2);
	gsl_matrix_memcpy(myself,&rhs());
	KILLME = true;
}

celerium::gsl::Matrix::Matrix(size_t n, const std::initializer_list<double> init){
	myself = gsl_matrix_alloc(n,n);
	KILLME = true;

	auto iter = init.begin();
	for(size_t i=0; i<n; ++i)
	  for(size_t j=0; j<n; ++j){
	    if(iter != init.end()){
	      gsl_matrix_set(myself,i,j,*iter);
	      iter++;
	    }
	    else gsl_matrix_set(myself,i,j,0.0);
	  }
}

namespace celerium{ namespace gsl{
std::ostream& operator<<(std::ostream& str, const celerium::gsl::Matrix& A){
	for(size_t i=0; i<A.rowNumber(); ++i){
	  str<<"| ";
	  for(size_t j=0; j<A.columnNumber(); ++j)
	    str<<A(i,j)<<" ";
	  str<<"|\n";
	}

	return str;
}
}}

celerium::gsl::Matrix::~Matrix(){
	if(KILLME) gsl_matrix_free(myself);
}

void celerium::gsl::Matrix::resize(size_t n, bool INITIALIZE){
	if(KILLME)	gsl_matrix_free(myself);
	if(INITIALIZE)  myself = gsl_matrix_calloc(n,n);
	else		myself = gsl_matrix_alloc(n,n);
	KILLME = true;
}

const gsl_matrix& celerium::gsl::Matrix::operator()() const{
	return *myself;
}

gsl_matrix* celerium::gsl::Matrix::operator()(){
	return myself; //unsafe!
}

double& celerium::gsl::Matrix::operator()(size_t i, size_t j){
	return *gsl_matrix_ptr(myself,i,j);
}

double celerium::gsl::Matrix::operator()(size_t i, size_t j) const{
	return *gsl_matrix_ptr(myself,i,j);
}

int celerium::gsl::Matrix::swap(celerium::gsl::Matrix& rhs){
	return gsl_matrix_swap(myself, rhs());
}

void celerium::gsl::Matrix::zero(){
	gsl_matrix_set_zero(myself);
}

void celerium::gsl::Matrix::one(){
	gsl_matrix_set_identity(myself);
}

size_t celerium::gsl::Matrix::rowNumber() const{
	return myself->size1;
}

size_t celerium::gsl::Matrix::columnNumber() const{
	return myself->size2;
}

celerium::gsl::Vector celerium::gsl::Matrix::row(size_t i){
	gsl_vector* tmp = gsl_vector_alloc(columnNumber());
	gsl_matrix_get_row(tmp, myself, i);

	return Vector(tmp);
}

celerium::gsl::Vector celerium::gsl::Matrix::column(size_t i){
	gsl_vector* tmp = gsl_vector_alloc(rowNumber());
	gsl_matrix_get_col(tmp, myself, i);

	return Vector(tmp);
}

void celerium::gsl::Matrix::invert(){
	if( rowNumber() != columnNumber())
	 throw std::logic_error("I'm unable to calculate Pfaffians! Please just work on Quadratic Matrices...");
	gsl_matrix* tmp = gsl_matrix_alloc(rowNumber(),rowNumber());
	gsl_permutation* permutation = gsl_permutation_alloc(rowNumber());
	int sign;

	gsl_linalg_LU_decomp(myself,permutation,&sign);
	gsl_linalg_LU_invert(myself,permutation,tmp);
	gsl_permutation_free (permutation);

	gsl_matrix_free(myself);
	myself = tmp;
}

int celerium::gsl::Matrix::transpose(){
	return gsl_matrix_transpose(myself);
}

/*
 * ARITHM. OPERATORS
 */

celerium::gsl::Matrix& celerium::gsl::Matrix::operator+=(const celerium::gsl::Matrix& rhs){
	if(rowNumber() != rhs.rowNumber() || columnNumber() != rhs.columnNumber())
	  throw std::logic_error("Oy! You can't add different-type matrices!");
	gsl_matrix_add(myself,&rhs());
	return *this;
}

celerium::gsl::Matrix celerium::gsl::Matrix::operator+(const celerium::gsl::Matrix& rhs) const{
	Matrix result(*this);
	return result+=rhs;
}

celerium::gsl::Matrix& celerium::gsl::Matrix::operator-=(const celerium::gsl::Matrix& rhs){
	if(rowNumber() != rhs.rowNumber() || columnNumber() != rhs.columnNumber())
	  throw std::logic_error("Oy! You can't substract different-type matrices!");
	gsl_matrix_sub(myself,&rhs());
	return *this;
}

celerium::gsl::Matrix celerium::gsl::Matrix::operator-(const celerium::gsl::Matrix& rhs) const{
	Matrix result(*this);
	return result-=rhs;
}

celerium::gsl::Matrix celerium::gsl::Matrix::operator*(const celerium::gsl::Matrix& rhs) const{
	if(rowNumber() != rhs.columnNumber())
	  throw std::logic_error("Oy! You can't multiply a ("
		 + std::to_string(rowNumber()) + "," 
		 + std::to_string(columnNumber()) + ") matrix by a (" 
		 + std::to_string(rhs.rowNumber()) + ","  
		 + std::to_string(rhs.columnNumber()) + ") matrix!");
	Matrix result(rhs.rowNumber(),columnNumber(),true);

	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,
			1.0,myself,&rhs(),
			0.0,result());

	return result;
}

celerium::gsl::Vector celerium::gsl::Matrix::operator*(const celerium::gsl::Vector& rhs) const{
	if(rowNumber() != rhs.size())
	  throw std::logic_error("Oy! You can't multiply a ("
		 + std::to_string(rowNumber()) + "," 
		 + std::to_string(columnNumber()) + ") matrix by a (" 
		 + std::to_string(rhs.size()) + ") vector!");
	Vector result(columnNumber(),true);

	gsl_blas_dgemv(CblasNoTrans, 1.0, myself, &rhs(),
			0.0, result());

	return result;
}

namespace celerium{ namespace gsl{
celerium::gsl::Vector operator*(const celerium::gsl::Vector& lhs, const celerium::gsl::Matrix& rhs){
	if(rhs.columnNumber() != lhs.size())
	  throw std::logic_error("Oy! You can't multiply a ("
		 + std::to_string(rhs.rowNumber()) + "," 
		 + std::to_string(rhs.columnNumber()) + ") matrix by a (" 
		 + std::to_string(lhs.size()) + ") vector!");
	celerium::gsl::Vector result(rhs.rowNumber(),true);

	gsl_blas_dgemv(CblasTrans, 1.0, &rhs(), &lhs(),
			0.0, result());

	return result;
}
}}

celerium::gsl::Matrix& celerium::gsl::Matrix::operator*=(const celerium::gsl::Matrix& rhs){
	if(rowNumber() != rhs.columnNumber())
	  throw std::logic_error("Oy! You can't multiply a ("
		 + std::to_string(rowNumber()) + "," 
		 + std::to_string(columnNumber()) + ") matrix by a (" 
		 + std::to_string(rhs.rowNumber()) + ","  
		 + std::to_string(rhs.columnNumber()) + ") matrix!");

	gsl_matrix* tmp = gsl_matrix_alloc(rowNumber(),columnNumber());
	gsl_matrix_memcpy(tmp,myself);

	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,
			1.0,tmp,&rhs(),
			0.0,myself);

	gsl_matrix_free(tmp);

	return *this;
}


celerium::gsl::Matrix& celerium::gsl::Matrix::operator=(const celerium::gsl::Matrix& rhs){
	if(KILLME) gsl_matrix_free(myself);
	myself = gsl_matrix_alloc(rhs().size1,rhs().size2);
	gsl_matrix_memcpy(myself,&rhs());
	KILLME = true;

	return *this;
}

celerium::gsl::Matrix& celerium::gsl::Matrix::operator=(const celerium::gsl::Vector& rhs){
	if(KILLME) gsl_matrix_free(myself);
	myself = gsl_matrix_calloc(rhs().size,rhs().size);
	
	for(size_t i=0; i<rhs().size; ++i)
		this->operator()(i,i) = rhs(i);

	KILLME = true;

	return *this;
}



/*
 * EIGENPROBLEM
 */

void celerium::gsl::Matrix::symmetricEigenProblem(celerium::gsl::Matrix& eigenvectors, celerium::gsl::Vector& eigenvalues){
	// makes sense only for symmetric matrices!
	if( rowNumber() != columnNumber())
	 throw std::logic_error("I'm unable to calculate Pfaffians! Please just work on Quadratic Matrices...");

	// resize output objects
	eigenvectors.resize(rowNumber());
	eigenvalues.resize(rowNumber());

	gsl_eigen_symmv_workspace *workspace = gsl_eigen_symmv_alloc(rowNumber());
	
	gsl_eigen_symmv(myself, eigenvalues(), eigenvectors(), workspace);

	eigenvectors.transpose(); // mathematical convention: vectors in columns!
	
	gsl_eigen_symmv_free(workspace);
}
