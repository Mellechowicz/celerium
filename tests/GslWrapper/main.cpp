#include "../../headers/gslmatrix.h"
#include <iostream>
#include <functional>

double pow4(const double& x){
	return pow(x,4);
}

struct Pow3{
	double operator()(const double& x){
		return pow(x,3);
	}
} pow3;

std::function<double(const double&)> pow2([](const double& x){return x*x;});

int main(){
	celerium::gsl::Matrix A(2,true);
	A(0,0) = 1.0;
	A(1,0) = 2.0;
	A(0,1) = 2.0;
	A(1,1) = 1.0;

	celerium::gsl::Matrix B(2,2,false);
	B.one();

	celerium::gsl::Matrix C(2,{2.0, -1.0,
			          -1.0});
	
	celerium::gsl::Matrix D = A;
	D.invert();

	std::cout<<"Algebra:"<<std::endl;
	std::cout<<"A:\n"<<A<<std::endl;
	std::cout<<"B:\n"<<B<<std::endl;
	std::cout<<"C:\n"<<C<<std::endl;
	std::cout<<"D:\n"<<D<<std::endl;
	std::cout<<"A+B:\n"<<A+B<<std::endl;
	std::cout<<"B+A:\n"<<B+A<<std::endl;
	std::cout<<"A-B:\n"<<A-B<<std::endl;
	std::cout<<"B-A:\n"<<B-A<<std::endl;
	std::cout<<"A*B:\n"<<A*B<<std::endl;
	std::cout<<"A*D:\n"<<A*D<<std::endl;


	std::cout<<"******************************************************************************"<<std::endl;
	std::cout<<"Eigenproblem:"<<std::endl;
	celerium::gsl::Matrix eigVecs;
	celerium::gsl::Vector eigVals(2);

	std::cout<<"For:\n"<<A;
	A.symmetricEigenProblem(eigVecs,eigVals);
	std::cout<<"Eigenvalues:\n"<<eigVals;
	std::cout<<"Eigenvectors:\n"<<eigVecs;

	celerium::gsl::Matrix eigVecsInv = eigVecs;
	eigVecsInv.invert();

	std::cout<<"After diagonlization:\n"<<eigVecs*A*eigVecsInv;


	std::cout<<"******************************************************************************"<<std::endl;
	std::cout<<"Applying functions for A:"<<std::endl<<A;

	std::cout<<"A^2 explicitly:\n"<<A*A;
	std::cout<<"A^2 optimally:\n"<<A.apply(pow2);
	std::cout<<"A^3 explicitly:\n"<<A*A*A;
	std::cout<<"A^3 optimally:\n"<<A.apply(pow3);
	std::cout<<"A^4 explicitly:\n"<<A*A*A*A;
	std::cout<<"A^4 optimally:\n"<<A.apply(pow4);
	std::cout<<"Trying sqrt(A):\n";
	try{
		std::cout<<A.apply([](double& x){return sqrt(x);});
	} catch(const std::exception& e){
		std::cerr<<"Unsuccessful: "<<e.what()<<std::endl;
	}

	std::cout<<"Trying 1/sqrt(A^2):\n";
	try{
		celerium::gsl::Matrix A2 = A*A;
		auto sqrtA2 = A2.apply([](double& x){return 1.0/sqrt(x);});
		std::cout<<sqrtA2;

		auto v0 = sqrtA2.column(0);
		auto v1 = sqrtA2.column(1);
		std::cout<<"v0: "<<v0;
		std::cout<<"v1: "<<v1;
	
		std::cout<<"v0 and v1 orthogonalized with respect to A2 = A*A"<<std::endl;

		std::cout<<"v0.v1 = "<<v0*v1<<std::endl;
		std::cout<<"v0.A2.v0 = "<<(v0*A)*(A*v0)<<std::endl;
		std::cout<<"v0.A2.v1 = "<<(v0*A)*(A*v1)<<std::endl;
		std::cout<<"v1.A2.v0 = "<<(v1*A)*(A*v0)<<std::endl;
		std::cout<<"v1.A2.v1 = "<<(v1*A)*(A*v1)<<std::endl;
		

	} catch(const std::exception& e){
		std::cerr<<"Unsuccessful: "<<e.what()<<std::endl;
	}

	return 0;
}
