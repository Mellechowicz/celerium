#include "../headers/gslmatrixcomplex.h"

celerium::gsl::VectorComplex::VectorComplex(){
	myself = nullptr;
	KILLME = false;
}

celerium::gsl::VectorComplex::VectorComplex(size_t n, bool INITIALIZE){
	KILLME = true;
	if(INITIALIZE)  myself = gsl_vector_complex_calloc(n);
	else		myself = gsl_vector_complex_alloc(n);
}

celerium::gsl::VectorComplex::VectorComplex(const VectorComplex& rhs){
	myself = gsl_vector_complex_alloc(rhs().size);
	gsl_vector_complex_memcpy(myself,&rhs());
	KILLME = true;
}

celerium::gsl::VectorComplex::VectorComplex(gsl_vector_complex* v){
	KILLME = true;
	myself = v;
}

celerium::gsl::VectorComplex::VectorComplex(const gsl_vector_complex* v){
	myself = gsl_vector_complex_alloc(v->size);
	gsl_vector_complex_memcpy(myself,v);
	KILLME = true;
}

namespace celerium{ namespace gsl{
std::ostream& operator<<(std::ostream& str, const celerium::gsl::VectorComplex& A){
	str<<"[ ";
	for(size_t i=0; i<A.size(); ++i)
	  str<<A(i)<<" ";
	str<<"]\n";

	return str;
}
}}

void celerium::gsl::VectorComplex::set(gsl_vector_complex* v){
	if(KILLME) gsl_vector_complex_free(myself);
	myself = v;
	KILLME = true;
}

celerium::gsl::VectorComplex::~VectorComplex(){
	if(KILLME) gsl_vector_complex_free(myself);
}

void celerium::gsl::VectorComplex::resize(size_t n, bool INITIALIZE){
	if(KILLME)	gsl_vector_complex_free(myself);
	if(INITIALIZE)  myself = gsl_vector_complex_calloc(n);
	else		myself = gsl_vector_complex_alloc(n);
	KILLME = true;
}

const gsl_vector_complex& celerium::gsl::VectorComplex::operator()() const{
	return *myself;
}

gsl_vector_complex* celerium::gsl::VectorComplex::operator()(){
	return myself; //unsafe!
}

int celerium::gsl::VectorComplex::swap(celerium::gsl::VectorComplex& rhs){
	return gsl_vector_complex_swap(myself, rhs());
}

std::complex<double> celerium::gsl::VectorComplex::operator()(size_t i) const{
  gsl_complex result = gsl_vector_complex_get(myself,i);
  return {GSL_REAL(result), GSL_IMAG(result)};
}
/*
double& celerium::gsl::VectorComplex::operator()(size_t i){
	return *gsl_vector_complex_ptr(myself,i);
}
*/
void celerium::gsl::VectorComplex::zero(){
  if (!KILLME) return; // Do nothing if matrix is not initialized;
	gsl_vector_complex_set_zero(myself);
}

int celerium::gsl::VectorComplex::one(size_t i){
	return gsl_vector_complex_set_basis(myself,i);
}

size_t celerium::gsl::VectorComplex::size() const{
	return myself->size;
}

/*
 * ARITHM. OPERATORS
 */

celerium::gsl::VectorComplex& celerium::gsl::VectorComplex::operator+=(const celerium::gsl::VectorComplex& rhs){
	if(size() != rhs.size())
	  throw std::logic_error("Oy! You can't add different-type vectors!");
	gsl_vector_complex_add(myself,&rhs());
	return *this;
}

celerium::gsl::VectorComplex celerium::gsl::VectorComplex::operator+(const celerium::gsl::VectorComplex& rhs) const{
	VectorComplex result(*this);
	return result+=rhs;
}

celerium::gsl::VectorComplex& celerium::gsl::VectorComplex::operator-=(const celerium::gsl::VectorComplex& rhs){
	if(size() != rhs.size())
	  throw std::logic_error("Oy! You can't substract different-type vectors!");
	gsl_vector_complex_sub(myself,&rhs());
	return *this;
}

celerium::gsl::VectorComplex celerium::gsl::VectorComplex::operator-(const celerium::gsl::VectorComplex& rhs) const{
	VectorComplex result(*this);
	return result-=rhs;
}

std::complex<double> celerium::gsl::VectorComplex::operator*(const celerium::gsl::VectorComplex& rhs) const{
	gsl_complex result;
	gsl_blas_zdotc(myself,&rhs(),&result);
	return {GSL_REAL(result), GSL_IMAG(result)};
}

celerium::gsl::MatrixComplex::MatrixComplex(){
	myself = nullptr;
	KILLME = false;
}

celerium::gsl::MatrixComplex::MatrixComplex(size_t n, bool INITIALIZE){
	  if(INITIALIZE) myself = gsl_matrix_complex_calloc(n,n);
	  else		 myself = gsl_matrix_complex_alloc(n,n);
	  KILLME = true;
}

celerium::gsl::MatrixComplex::MatrixComplex(size_t n, size_t m, bool INITIALIZE){
	if(INITIALIZE)  myself = gsl_matrix_complex_calloc(n,m);
	else		myself = gsl_matrix_complex_alloc(n,m);
	KILLME = true;
}

celerium::gsl::MatrixComplex::MatrixComplex(const celerium::gsl::MatrixComplex& rhs){
	myself = gsl_matrix_complex_alloc(rhs().size1,rhs().size2);
	gsl_matrix_complex_memcpy(myself,&rhs());
	KILLME = true;
}

celerium::gsl::MatrixComplex::MatrixComplex(size_t n, const std::initializer_list<std::complex<double>> init){
	myself = gsl_matrix_complex_alloc(n,n);
	KILLME = true;

	auto iter = init.begin();
	for(size_t i=0; i<n; ++i)
	  for(size_t j=0; j<n; ++j){
	    if(iter != init.end()){
	      gsl_matrix_complex_set(myself,i,j,{(*iter).real(), (*iter).imag()});
	      iter++;
	    }
	    else gsl_matrix_complex_set(myself,i,j,{0.0, 0.0});
	  }
}

namespace celerium{ namespace gsl{
std::ostream& operator<<(std::ostream& str, const celerium::gsl::MatrixComplex& A){
	for(size_t i=0; i<A.rowNumber(); ++i){
	  str<<"| ";
	  for(size_t j=0; j<A.columnNumber(); ++j)
	    str<<A(i,j)<<" ";
	  str<<"|\n";
	}

	return str;
}
}}

celerium::gsl::MatrixComplex::~MatrixComplex(){
	if(KILLME) gsl_matrix_complex_free(myself);
}

void celerium::gsl::MatrixComplex::resize(size_t n, bool INITIALIZE){
	if(KILLME)	gsl_matrix_complex_free(myself);
	if(INITIALIZE)  myself = gsl_matrix_complex_calloc(n,n);
	else		myself = gsl_matrix_complex_alloc(n,n);
	KILLME = true;
}

const gsl_matrix_complex& celerium::gsl::MatrixComplex::operator()() const{
	return *myself;
}

gsl_matrix_complex* celerium::gsl::MatrixComplex::operator()(){
	return myself; //unsafe!
}

/*
double& celerium::gsl::MatrixComplex::operator()(size_t i, size_t j){
	return *gsl_matrix_complex_ptr(myself,i,j);
}
*/

double& celerium::gsl::MatrixComplex::real(size_t i, size_t j) {
  return GSL_REAL(*gsl_matrix_complex_ptr(myself,i,j));
}
double& celerium::gsl::MatrixComplex::imag(size_t i, size_t j) {
  return GSL_IMAG(*gsl_matrix_complex_ptr(myself,i,j));
}


std::complex<double> celerium::gsl::MatrixComplex::operator()(size_t i, size_t j) const{
        gsl_complex &result = *gsl_matrix_complex_ptr(myself,i,j);
	return {GSL_REAL(result), GSL_IMAG(result)};
}

int celerium::gsl::MatrixComplex::swap(celerium::gsl::MatrixComplex& rhs){
	return gsl_matrix_complex_swap(myself, rhs());
}

void celerium::gsl::MatrixComplex::zero(){
    if (!KILLME) return; // Do nothing if matrix is not initialized;
	gsl_matrix_complex_set_zero(myself);
}

void celerium::gsl::MatrixComplex::one(){
	gsl_matrix_complex_set_identity(myself);
}

size_t celerium::gsl::MatrixComplex::rowNumber() const{
	return myself->size1;
}

size_t celerium::gsl::MatrixComplex::columnNumber() const{
	return myself->size2;
}

celerium::gsl::VectorComplex celerium::gsl::MatrixComplex::row(size_t i){
	gsl_vector_complex* tmp = gsl_vector_complex_alloc(columnNumber());
	gsl_matrix_complex_get_row(tmp, myself, i);

	return VectorComplex(tmp);
}

celerium::gsl::VectorComplex celerium::gsl::MatrixComplex::column(size_t i){
	gsl_vector_complex* tmp = gsl_vector_complex_alloc(rowNumber());
	gsl_matrix_complex_get_col(tmp, myself, i);

	return VectorComplex(tmp);
}

void celerium::gsl::MatrixComplex::invert(){
	if( rowNumber() != columnNumber())
	 throw std::logic_error("I'm unable to calculate Pfaffians! Please just work on Quadratic Matrices...");
	gsl_matrix_complex* tmp = gsl_matrix_complex_alloc(rowNumber(),rowNumber());
	gsl_permutation* permutation = gsl_permutation_alloc(rowNumber());
	int sign;

	gsl_linalg_complex_LU_decomp(myself,permutation,&sign);
	gsl_linalg_complex_LU_invert(myself,permutation,tmp);
	gsl_permutation_free (permutation);

	gsl_matrix_complex_free(myself);
	myself = tmp;
}

int celerium::gsl::MatrixComplex::transpose(){
	return gsl_matrix_complex_transpose(myself);
}

/*
 * ARITHM. OPERATORS
 */

celerium::gsl::MatrixComplex& celerium::gsl::MatrixComplex::operator+=(const celerium::gsl::MatrixComplex& rhs){
	if(rowNumber() != rhs.rowNumber() || columnNumber() != rhs.columnNumber())
	  throw std::logic_error("Oy! You can't add different-type matrices!");
	gsl_matrix_complex_add(myself,&rhs());
	return *this;
}

celerium::gsl::MatrixComplex celerium::gsl::MatrixComplex::operator+(const celerium::gsl::MatrixComplex& rhs) const{
	MatrixComplex result(*this);
	return result+=rhs;
}

celerium::gsl::MatrixComplex& celerium::gsl::MatrixComplex::operator-=(const celerium::gsl::MatrixComplex& rhs){
	if(rowNumber() != rhs.rowNumber() || columnNumber() != rhs.columnNumber())
	  throw std::logic_error("Oy! You can't substract different-type matrices!");
	gsl_matrix_complex_sub(myself,&rhs());
	return *this;
}

celerium::gsl::MatrixComplex celerium::gsl::MatrixComplex::operator-(const celerium::gsl::MatrixComplex& rhs) const{
	MatrixComplex result(*this);
	return result-=rhs;
}

celerium::gsl::MatrixComplex celerium::gsl::MatrixComplex::operator*(const celerium::gsl::MatrixComplex& rhs) const{
	if(rowNumber() != rhs.columnNumber())
	  throw std::logic_error("Oy! You can't multiply a ("
		 + std::to_string(rowNumber()) + "," 
		 + std::to_string(columnNumber()) + ") matrix by a (" 
		 + std::to_string(rhs.rowNumber()) + ","  
		 + std::to_string(rhs.columnNumber()) + ") matrix!");
	MatrixComplex result(rhs.rowNumber(),columnNumber(),true);

	gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,
                       {1.0,0.0},myself,&rhs(),
                       {0.0,0.0},result());

	return result;
}

celerium::gsl::VectorComplex celerium::gsl::MatrixComplex::operator*(const celerium::gsl::VectorComplex& rhs) const{
	if(rowNumber() != rhs.size())
	  throw std::logic_error("Oy! You can't multiply a ("
		 + std::to_string(rowNumber()) + "," 
		 + std::to_string(columnNumber()) + ") matrix by a (" 
		 + std::to_string(rhs.size()) + ") vector!");
	VectorComplex result(columnNumber(),true);

	gsl_blas_zgemv(CblasNoTrans, {1.0, 0.0}, myself, &rhs(),
                       {0.0,0.0}, result());

	return result;
}

namespace celerium{ namespace gsl{
celerium::gsl::VectorComplex operator*(const celerium::gsl::VectorComplex& lhs, const celerium::gsl::MatrixComplex& rhs){
	if(rhs.columnNumber() != lhs.size())
	  throw std::logic_error("Oy! You can't multiply a ("
		 + std::to_string(rhs.rowNumber()) + "," 
		 + std::to_string(rhs.columnNumber()) + ") matrix by a (" 
		 + std::to_string(lhs.size()) + ") vector!");
	celerium::gsl::VectorComplex result(rhs.rowNumber(),true);

	gsl_blas_zgemv(CblasTrans, {1.0, 0.0}, &rhs(), &lhs(),
                       {0.0,0.0}, result());

	return result;
}
}}

celerium::gsl::MatrixComplex& celerium::gsl::MatrixComplex::operator*=(const celerium::gsl::MatrixComplex& rhs){
	if(rowNumber() != rhs.columnNumber())
	  throw std::logic_error("Oy! You can't multiply a ("
		 + std::to_string(rowNumber()) + "," 
		 + std::to_string(columnNumber()) + ") matrix by a (" 
		 + std::to_string(rhs.rowNumber()) + ","  
		 + std::to_string(rhs.columnNumber()) + ") matrix!");

	gsl_matrix_complex* tmp = gsl_matrix_complex_alloc(rowNumber(),columnNumber());
	gsl_matrix_complex_memcpy(tmp,myself);

	gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,
                       {1.0,0.0},tmp,&rhs(),
                       {0.0,0.0},myself);

	gsl_matrix_complex_free(tmp);

	return *this;
}


celerium::gsl::MatrixComplex& celerium::gsl::MatrixComplex::operator=(const celerium::gsl::MatrixComplex& rhs){
	if(KILLME) gsl_matrix_complex_free(myself);
	myself = gsl_matrix_complex_alloc(rhs().size1,rhs().size2);
	gsl_matrix_complex_memcpy(myself,&rhs());
	KILLME = true;

	return *this;
}

celerium::gsl::MatrixComplex& celerium::gsl::MatrixComplex::operator=(const celerium::gsl::VectorComplex& rhs){
	if(KILLME) gsl_matrix_complex_free(myself);
	myself = gsl_matrix_complex_calloc(rhs().size,rhs().size);
	
	for(size_t i=0; i<rhs().size; ++i)
		this->operator()(i,i) = rhs(i);

	KILLME = true;

	return *this;
}

celerium::gsl::MatrixComplex& celerium::gsl::MatrixComplex::operator*(std::complex<double> scalar) {
  gsl_matrix_complex_scale(myself,{scalar.real(),scalar.imag()});
  return *this;
}



/*
 * EIGENPROBLEM
 */

void celerium::gsl::MatrixComplex::symmetricEigenProblem(celerium::gsl::MatrixComplex& eigenvectors, celerium::gsl::Vector& eigenvalues){
	// makes sense only for symmetric matrices!
	if( rowNumber() != columnNumber())
	 throw std::logic_error("I'm unable to calculate Pfaffians! Please just work on Quadratic Matrices...");

	// resize output objects
	eigenvectors.resize(rowNumber());
	eigenvalues.resize(rowNumber());

	gsl_eigen_hermv_workspace *workspace = gsl_eigen_hermv_alloc(rowNumber());
	
	gsl_eigen_hermv(myself, eigenvalues(), eigenvectors(), workspace);

	eigenvectors.transpose(); // mathematical convention: vectors in columns!
	
	gsl_eigen_hermv_free(workspace);
}
